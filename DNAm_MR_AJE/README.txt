README
Documentation for the paper "DNA methylation as a possible causal mechanism linking childhood adversity and health: Results from two-sample mendelian randomization study." published in American Journal of Epidemiology, in the Special Collection for Mental Health
Authors: 	Isabel K. Schuurmans; Erin C. Dunn; Alexandre A. Lussier
Last updated:	March 8, 2024

In this folder, you can find the following information regarding our manuscript. 

1. Manuscripts
	- Final version of the manuscript and supplementary appendices, before journal editing 

2. Scripts 
R scripts used to analyze each adversity for each of the main analyses using the Structured Life Course Modeling Approach (SLCMA)
        - runMRfunction_22024-03-08.R - the function needed to run MR analyses
	- 1. prepareMR_2024-03-08.R - a script to prepare summary statistics for MR
	- 2. runningMR_2024-03-08.R - a script to run the mendelian randomization; generating raw results
	- 3. analyzingMR_2024-03-08.R - a script to analyze raw MR results, and to produce sensitivity checks
        - 4. createTablesMR_2024-03-08.R - a script to generate tables and figures

Summary statistics used: The GWAS summary statistics used for the analysis were downloaded from publicly available sources. 
- Summary statistics for DNAm loci were downloaded from http://godmc.org.uk/. 
- Summary statistics for ADHD, anorexia nervosa, anxiety disorders, autism spectrum disorder, bipolar disorder, cannabis use disorder, obsessive-compulsive disorder, opioid exposed, post-traumatic stress disorder (PTSD), schizophrenia, and Tourette syndrome can be downloaded from https://pgc.unc.edu/. 
- Summary statistics for depression can be downloaded from https://ipsych.dk/en/research/downloads. 
- Summary statistics from suicide attempt can be downloaded from https://tinyurl.com/ISGC2021. 
- Summary statistics for stroke, asthma, and COPD were downloaded from https://www.globalbiobankmeta.org/resources. 
- Summary statistics for chronic kidney disease can be downloaded from http://ckdgen.imbi.uni-freiburg.de/. 
- Summary statistics for coronary artery disease were downloaded from https://cvd.hugeamp.org/downloads.html#summary. 
- Summary statistics for obesity can be downloaded from https://gwas.mrcieu.ac.uk/datasets/ukb-b-15541. 
- Summary statistics for type 2 diabetes can be downloaded from http://cnsgenomics.com/data.html. 
- Summary statistics for smoking initiation can be downloaded from https://gwas.mrcieu.ac.uk/datasets/ieu-b-4877. 
- Summary statistics for physical in activity were downloaded from https://www.ebi.ac.uk/gwas/publications/36071172. 

Please reach out to alussier[at]mgh.harvard.edu or edunn2[at]mgh.harvard.edu further details of these analyses. 