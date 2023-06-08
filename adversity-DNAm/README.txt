README
Documentation for the paper "Updates to data versions and analytic methods influence the reproducibility of results from epigenome-wide association studies"
Author: 	Alexandre Lussier, PhD
Last updated:	June 8, 2023

In this folder, you can find the following information regarding our manuscript. 

1. Manuscripts
	- Preprint of the manuscript by Lussier and colleagues described here, and supplemental materials (In press at Lancet Child and Adolescent Health).
	- Manuscript by Zhu and colleagues, which describes the application of the SLCMA to epigenome-wide settings. [Link to paper](https://static1.squarespace.com/static/62d6cd5bed9f7007dc679802/t/62faa96977092f5ee08c2e9f/1660594538637/2021-06+Zhu+%28American+Journal+of+Epidemiology%29++A+Structured+Approach+to+Evaluating+Life-Course+Hypotheses-+Moving+Beyond+Analyses+of+Exposed+Versus+Unexposed+in+the+-Omics+Context.pdf)
	- Manuscript by Dunn and colleagues, which describes the application of the SLCMA to analyze the time-dependent relationship between childhood adversity and DNA methylation. [Link to paper](https://static1.squarespace.com/static/62d6cd5bed9f7007dc679802/t/62faadbac25ee6601e977719/1660595656020/2019-05+Dunn+%28Biol+Psych%29+-+Sensitive+Periods+for+the+Effect+of+Childhood+Adversity+on+DNA+Methylation-+Results+From+a+Prospective%2C+Longitudinal+Study.pdf)
	- Manuscript by Lussier and colleagues, which describes the findings between childhood adversity and DNA methylation at age 7. [Link to paper on results](https://static1.squarespace.com/static/62d6cd5bed9f7007dc679802/t/62faa70e0db3851fbf8ccae6/1660593937694/2022-02+Lussier+%28Epigenetics%29+Updates+to+data+versions+and+analytic+methods+influence+the+reproducibility+of+results+from+epigenome-wide+association+studies.pdf) and [link to paper on biological interpretation](https://static1.squarespace.com/static/62d6cd5bed9f7007dc679802/t/6307b1c699e7f30bb1c6e872/1661448646838/Lussier+2022+%28Bio+Psych+GOS%29+-+Sensitive+periods+for+the+effect+of+childhood+adversity+on+DNA+methylaton-+Updated+results+from+prospective%2C+longitudinal+study.pdf)


2. SLCMA scripts 
R scripts used to analyze each adversity for each of the main analyses using the Structured Life Course Modeling Approach (SLCMA)
	- 1. LARS-noimpute-function-20190516 - SLCMA function used for these analyses
	- 2. runSLCMA.adversity.noever.20200513.R - primary script, used for abuse, mompsych, oneadult, and nbhqual
	- 3. 2021-04-24_r_faminstability_slcma_noever.R - same as primary, but using r_faminstability 
	- 4. 2021-06-01_fscore_slcma_noever.R - same as primary, but using revised Fscore variable
	- 5. 2022-07-13_parcruelty_revised_SLCMA_noever_SI.R - same as primary, but using revised parcruelty variable

3. Analysis scripts
Scripts to analyze the results from the SLCMA and generate the figures
	- 1. Results_structure_2022-08-18.Rmd - concatenate results from the SLCMA scripts 
	- 2. Adversity_summary_2022-08-18.Rmd - summary of childhood adversity measures (Fig 1)
	- 3. Top_hits_biology_2022-08-18.Rmd - analysis of main age 15 findings 
	- 4. Top_hits_stability_2022-12-05.Rmd - stability of DNAm changes from age 7 and 15
	- 5. Top_hits_bootstrap_2022-08-24.Rmd - internal validation using non-parametric bootstrap
	- 6. Top_hits_confounders_2022-12-03.Rmd - sensitivity analyses of confounders and mediators
	- 7. Top_hits_MAD_models_2022-08-26.Rmd - sensitivity analyses adjusting for exposures to other adversities
	- 8. Top_hits_trajectories_2022-09-02.Rmd - analysis of the types of DNAm trajectories
	- 9. Additional_analyses_2022-12-09.Rmd - from revisions (QQ plots; stability of DNAm without adversity; EWAS catalog lookup; missMethyl gene ontology)



Please reach out to alussier[at]mgh.harvard.edu or edunn2[at]mgh.harvard.edu for summary statistics or further details of these analyses. 

