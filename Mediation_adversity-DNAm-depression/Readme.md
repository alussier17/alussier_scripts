README
Documentation for the paper "DNA methylation mediates the link between adversity and depressive symptoms" published in Nature Mental Health (2024)
Author: 	Alexandre Lussier, PhD
Last updated:	May 5, 2024

In this folder, you can find the following information regarding our manuscript. 

1. Manuscripts
 - Primary manuscript
 - Supplemental materials
 - Research briefing summarizing the main findings
 - Editorial article on the manuscript

2. Primary analyses
 - TE_SLCMA_2021-07-29.Rmd -	SLCMA to identify life course hypothesis for each of 7 adversities.
 - create_regressed_mvalues_2021-03-25.R	- Script to prepare DNAm data for mediation analysis. Regressed all covariates on DNAm values (x7).
 - Schaid_med_fn_part1_2021-03-25.R	- Runs Schaid method through the creation of fit.grid
 - Schaid_med_fn_part2_2021-03-26.R	- Uses grids and inputs from part1 to derive effect estimates and standard errors for each mediator using regmed and lavaan packages.
 - Meds_pval+CI_creation_2021-06-06.R	- Script to derive pvalues and CIs using MCM

3. Replication
 - FFCWS_site_replication_2022-07-11B.Rmd - base script for replication in FFCWS
 - GenR_site_replication_2022-09-07.Rmd - base script for replication in GenR
 - Mediation_replication_newFig_2024-08-27.Rmd - script for primary replication figure
 - MediationResults_2023-06-19.Rmd - compilation of mediation replication results

4. Revisions
 - mediate_function.R - function for mediation analyses
 - mediation_revisions_2024-03-25.Rmd - revisions script 
