README of mediation scripts
Filename				Purpose

2. TE_SLCMA_2021-07-29.Rmd	SLCMA to identify life course hypothesis for each of 7 adversities. 
	-Mannan, this step will probably be omitted in your own analysis, but included for reference. 

3. create_regressed_mvalues_2021-03-25.R	Script to prepare DNAm data for mediation analysis. Regressed all covariates on DNAm values (x7).

4. Schaid_med_fn_part1_2021-03-25.R	Runs Schaid method through the creation of fit.grid

5. Schaid_med_fn_part2_2021-03-26.R	Uses grids and inputs from part1 to derive effect estimates and standard errors for each mediator using regmed and lavaan packages.

6. Meds_pval+CI_creation_2021-06-06.R	Script to derive pvalues and CIs using MCM

