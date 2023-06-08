README
Documentation for the paper "Updates to data versions and analytic methods influence the reproducibility of results from epigenome-wide association studies"
Author: 	Alexandre Lussier, PhD
Last updated:	May 10, 2023

In this folder, you can find the following information regarding our manuscript. 

0. Manuscript
	- Main article and supplemental material 
	- Correspondence article describing similarities and differences with prior findings by Dunn et al.

1. Related articles 
	- Manuscript by Zhu and colleagues, which describes the application of the SLCMA to epigenome-wide settings.
	- Manuscript by Dunn and colleagues, which is the precursor to this paer.
	
	- Manuscript by Dunn and colleagues, which describes the application of the SLCMA to analysis the time-dependent relationship between childhood adversity and DNA methylation.
	- Supplemental

2. SLCMA scripts 
	- R scripts used to analyze each adversity for each of the main analyses
	- 01-oldDNAm-origCovars-covTest - old DNAm data, standard covariates, covariance test
	- 02-newDNAm-origCovars-covTest - new DNAm data, standard covariates, covariance test
	- 03-newDNAm-origCovars-sI - new DNAm data, standard covariates, selective inference
	- 04-newDNAm-educFWL-sI - new DNAm data, revised covariates, FWL adjustment, selective inference


3. Scripts to analyze the resulting data and generate the figures
	- 1. DD_overall_differences_2021-09-21.Rmd - general differences between data versions
	- 2. DD_EWAS_2021-08-20.Rmd - epigenome-wide analysis of childhood adversity or prenatal smoking
	- 3. DD_SLCMA_process_2021-09-21.Rmd - concatenate results from the SLCMA scripts
	- 4. DD_figures_2021-08-26.Rmd - final analysis and figures for publication
