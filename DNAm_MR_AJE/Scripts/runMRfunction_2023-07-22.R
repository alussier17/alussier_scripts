#' This function estimates the causal relationship between DNA methylation (DNAm) and health outcomes using single nucleotide polymorphisms (SNPs) as instrumental variables (IVs), emplying the TwoSampleMR package.
#'
#' @param cpgs Vector specifying all CpGs of interest. This must be specified.
#' @param age When using 'mQTLdb' data, specify the age for the data usage. Use 'Childhood' for age 7, and 'Adolescent' for age 15. This is required if 'mQTLdb' is used.
#' @param database Specify which mQTL database to use: 'mQTLdb', or 'goDMC'. This must be specified.
#' @param clumpDMC Logical, indicating whether to perform clumping on DMC. Default is TRUE.
#' @param psychopathologies Vector of outcomes to investigate. Defaults to an empty string. Ensure the list matches the names of GWAS summary statistics files intended for use. Defaults to "all". 
#' @param clumpMR Logical, indicating whether to exclude palindromic SNPs. Defaults to FALSE.
#' @param cis.pvalue P-value threshold for identifying SNPs associated with DNAm in cis (<1 Mb from loci). Defaults to 1e-8.
#' @param trans.pvalue P-value threshold for identifying SNPs associated with DNAm in trans (>1 Mb from loci). Defaults to 1e-14.
#' @param clump.r2 Threshold for excluding SNPs with high linkage disequilibrium (R2 > 0.01). Defaults to 0.01.
#' @param transQTLs Logical, indicating whether to include transQTLs. Defaults to FALSE.
#' @param verbose Logical, indicating whether to print detailed progress information. Defaults to TRUE.
#' @param multiVariate Logical, indicating support for multivariate analysis. Currently not supported. Defaults to FALSE.
#' @param largeFiles Logical, indicating whether large files are to be used. Defaults to TRUE.
#' @param skipOutcome Vector of outcomes to be skipped. Defaults to an empty string.
#' @param name Character string specifying the output naming convention. Defaults to an empty string.
#' @param path_to_stats Character string specifying the path to the GWAS summary statistics. Defaults to a predefined path in the documentation.
#' @param stat_summary Character string specifying the RData file summarizing all outcomes, prepared in a specified script. Defaults to "main_outcomes_2023-04-19.Rdata". The name of the object loaded in MUST be psych.
#' @param outputdir Character string specifying the directory for saving output data. Defaults to a predefined path in the documentation.
#'
#' @details The MR analysis relies on three critical assumptions:
#' 1. SNPs selected as IVs are strongly associated with DNAm.
#' 2. SNPs as IVs are not associated with confounders of the DNAm-health outcome association.
#' 3. SNPs affect health outcomes only through DNAm.
#' Violations of these assumptions may lead to biased causal estimates.
#'
#' @return The function returns an object containing the causal estimates of DNAm on health outcomes, along with various diagnostics related to the MR analysis.
#'
#' @examples
#' # Example usage of the runMR function:
#' results <- runMR(cpgs = "cg00000029", age = "Childhood", database = "mQTLdb",
#'                  psychopathologies = c("Schizophrenia", "Bipolar Disorder"),
#'                  path_to_stats = "~/path/to/summstats/",
#'                  stat_summary = "psych_stat_summary_2023-03-03.Rdata",
#'                  outputdir = "~/path/to/results/")
#'
#' @note It's important to ensure that all preconditions and assumptions of MR analysis are met before interpreting the results.


runMR <- function(cpgs = "NA", age = NA, database = NA, clumpDMC = T,
                  psychopathologies = "all", clumpMR=F,
                  cis.pvalue = 1e-8, trans.pvalue = 1e-14,
                  clump.r2 = 0.01, transQTLs = F, 
                  verbose=T, multiVariate =F, 
                  largeFiles = T, skipOutcome ="",
                  name = "",
                  path_to_stats = "Analysis/Data/Summstats/",
                  stat_summary = "main_outcomes_2023-04-19.Rdata",
                  outputdir = "Analysis/Results/"
){
  #error messages
  if(cpgs[1] == "NA")
  {stop("Please include a list of CpGs to analyze")}
  if(!database %in% c("mQTLdb", "goDMC"))
  {stop("Please specific the mQTL database to use from 'mQTLdb', or 'goDMC'")}
  if(is.na(age) & database %in% c("mQTLdb"))
  {stop("Please specify age at DNAm collection to analyze mQTLdb data (7, 'Childhood', 15, 'Adolescent')")}
  
  #runtime
  t <- Sys.time()
  print(paste0("Start time: ", t))
  
  #loading in list of summary stats
  x = load(paste0(path_to_stats, stat_summary))
  psych = get(x); rm(x)
  
  if(psychopathologies[1] =="")
  {stop(cat("Please select one or more psychopathologies from the following list (MR outcome):
            all, ", paste(psych$psych, collapse=", ")),"\n")}
  
  if(!require(TwoSampleMR)){remotes::install_github("MRCIEU/TwoSampleMR@0.4.26")} 
  library(TwoSampleMR)
  library(dplyr)
  library(data.table)
  
  #minor formatting
  cpgs = as.character(cpgs)
  #making working directory
  cat("Creating working directory \n")
  ndirs <- length(grep(paste0(name,"_", Sys.Date()), list.dirs(outputdir)))+1
  newDir <- paste0(Sys.Date(),toupper(letters)[ndirs])
  dir.create(paste0(outputdir, paste0(name,"_"), newDir))
  setwd(paste0(outputdir, paste0(name,"_"), newDir))
  rm(ndirs)
  
  details <- rbind(cpgs = paste(cpgs,collapse = ", "), 
                   age = age, database = database, clumpDMC = clumpDMC,
                   psychopathologies = psychopathologies, 
                   cis.pvalue = cis.pvalue, 
                   transQTLs = transQTLs, 
                   trans.pvalue = trans.pvalue,
                   clumpMR = clumpMR, clump.r2 = clump.r2, 
                   multiVariate = multiVariate, largeFiles = largeFiles,
                   name = name, stat_summary = stat_summary,
                   path_to_stats = path_to_stats, outputdir=outputdir)
  write.table(details, file =paste0("runMR_call.txt"), quote=F, 
              row.names=T, col.names = F, sep="\t")
  rm(details)
  
  #making record of output
  zz <- file("output.log", open = "wt")
  sink(zz ,type = "output", split=T)
  sink(zz, split=T)
  
  print(paste0("Start time: ", t))
  
  cat("Step 1 - reading in mQTL data \n")
  if(database=="mQTLdb"){
    if(age %in% c(7, "Childhood")){
      cat("Obtaining mQTLs from mQTLdb - childhood DNAm \n")
      db <- fread(paste0(path_to_stats,"mQTLs/mQTLdb.F7.ALL.M.tab"))
      db <- db[db$gene %in% cpgs,]
      
      cat("Identifying cis and trans mQTLs \n")
      cgmap <- fread(paste0(path_to_stats,"mQTLs/ariesmqtlcpgs.csv"))
      cgmap <- cgmap[match(db$gene, cgmap$TargetID),]
      snpmap <- fread(paste0(path_to_stats,"mQTLs/ariesmqtlsnps.bim"))
      snpmap <- snpmap[match(db$SNP,snpmap$V2), ]
      head(snpmap)
      db <- cbind(db,
                  SNP_chr = snpmap$V1, 
                  SNP_pos = snpmap$V4,
                  cpg_chr = cgmap$CHR37,
                  cpg_pos = cgmap$COORDINATE_37,
                  A1 = snpmap$V5,
                  A2 = snpmap$V6)
      db$trans <- ifelse(db$SNP_chr != db$cpg_chr, "Y",
                         ifelse(abs(db$SNP_pos - db$cpg_pos) >1e6, "Y","N"))
      rm(cgmap, snpmap)
      
      exposureSNPs <- db[db$trans=="N" & db$`p-value`<= cis.pvalue,]
      if(transQTLs ==T){
        exposureSNPs2 <- db[db$trans=="Y" & db$`p-value`<= trans.pvalue,]
        exposureSNPs <- rbind(exposureSNPs, exposureSNPs2)
        rm(exposureSNPs2)
      }
      rm(db)
    }
    
    if(age %in% c(15, "Adolescent")){
      cat("Obtaining mQTLs from mQTLdb - adolescent DNAm \n")
      db <- data.table::fread(paste0(path_to_stats,"mQTLs/mQTLdb.15up.ALL.M.tab"))
      db <- db[db$gene %in% cpgs,]
      
      cat("Identifying cis and trans mQTLs \n")
      cgmap <- data.table::fread(paste0(path_to_stats,"mQTLs/ariesmqtlcpgs.csv"))
      cgmap <- cgmap[match(db$gene, cgmap$TargetID),]
      snpmap <- data.table::fread(paste0(path_to_stats,"mQTLs/ariesmqtlsnps.bim"))
      snpmap <- snpmap[match(db$SNP,snpmap$V2), ]
      
      db <- cbind(db,
                  SNP_chr = snpmap$V1, 
                  SNP_pos = snpmap$V4,
                  cpg_chr = cgmap$CHR37,
                  cpg_pos = cgmap$COORDINATE_37)
      db$trans <- ifelse(db$SNP_chr != db$cpg_chr, "Y",
                         ifelse(abs(db$SNP_pos - db$cpg_pos) >1e6, "Y","N"))
      rm(cgmap, snpmap)
      
      exposureSNPs <- db[db$trans=="N" & db$`p-value`<= cis.pvalue,]
      if(transQTLs ==T){
        exposureSNPs2 <- db[db$trans=="Y" & db$`p-value`<= trans.pvalue,]
        exposureSNPs <- rbind(exposureSNPs, exposureSNPs2)
        rm(exposureSNPs)
      }
      rm(db)
    }
    
    colnames(exposureSNPs)[colnames(exposureSNPs) =="gene"] <-"cpg"
    colnames(exposureSNPs)[colnames(exposureSNPs) =="p-value"] <-"p.value"
    exposureSNPs$SE <- exposureSNPs$beta/exposureSNPs$`t-stat`
    
  }
  
  if(database=="goDMC"){
    suppressMessages(library(jsonlite))
    suppressMessages(library(httr))
    suppressMessages(library(httr))
    if(clumpDMC ==T){
      cat("Obtaining mQTLs from goDMC - clumped \n")
      query <- list(
        cpgs = as.character(cpgs),
        pval = cis.pvalue,
        cistrans = "cis",
        clumped = 1
      )
      res <- POST("http://api.godmc.org.uk/v0.1/query", body = query, encode = "json")
      exposureSNPs <- content(res) %>% lapply(., as_data_frame) %>% bind_rows
      
      if(transQTLs ==T){
        cat("Including trans-mQTLs \n")
        query <- list(
          cpgs = as.character(cpgs),
          pval = trans.pvalue,
          cistrans = "trans",
          clumped = 1
        )
        res <- POST("http://api.godmc.org.uk/v0.1/query", body = query, encode = "json")
        exposureSNPs2 <- content(res) %>% lapply(., as_data_frame) %>% bind_rows
        exposureSNPs <- rbind(exposureSNPs, exposureSNPs2)
        rm(exposureSNPs2)
      }
      rm(query, res)
    }
    
    if(clumpDMC ==F){
      cat("Obtaining mQTLs from goDMC - unclumped \n")
      query <- list(
        cpgs = as.character(cpgs),
        pval = cis.pvalue,
        cistrans = "cis",
        clumped = 0
      )
      res <- POST("http://api.godmc.org.uk/v0.1/query", body = query, encode = "json")
      exposureSNPs <- content(res) %>% lapply(., as_data_frame) %>% bind_rows
      
      if(transQTLs ==T){
        cat("Including trans-mQTLs \n")
        query <- list(
          cpgs = as.character(cpgs),
          pval = trans.pvalue,
          cistrans = "trans",
          clumped = 0
        )
        res <- POST("http://api.godmc.org.uk/v0.1/query", body = query, encode = "json")
        exposureSNPs2 <- content(res) %>% lapply(., as_data_frame) %>% bind_rows
        exposureSNPs <- rbind(exposureSNPs, exposureSNPs2)
        rm(exposureSNPs2)
      }
      rm(query, res)
    }
    colnames(exposureSNPs)[1:2] <- toupper(colnames(exposureSNPs)[1:2])
    colnames(exposureSNPs)[colnames(exposureSNPs)=="rsid"] <- "SNP"
    colnames(exposureSNPs)[colnames(exposureSNPs)=="beta_are_a1"] <- "beta"
    colnames(exposureSNPs)[colnames(exposureSNPs)=="se_are"] <- "SE"
    colnames(exposureSNPs)[colnames(exposureSNPs)=="pval_are"] <- "p.value"
  }
  
  cat(paste0(nrow(exposureSNPs), " SNPs found for ", length(unique(exposureSNPs$cpg)), " CpGs \n"))
  save(exposureSNPs, file=paste0("exposureSNPs_", database, "_", newDir,".Rdata" ))
  write.csv(exposureSNPs, file = paste0("exposureSNPs_", database, "_", newDir,".csv"), quote=F)
  
  cat("\nStep 2 - MR of DNAm -> phenotype \n")
  dnam <- read_exposure_data(filename = paste0("exposureSNPs_", database, "_", newDir,".csv"),
                             sep=",",
                             snp_col = "SNP",
                             beta_col = "beta",
                             se_col = "SE",
                             effect_allele_col = "A1",
                             other_allele_col = "A2",
                             pval_col = "p.value",
                             gene_col="cpg")
  dnam$exposure <- dnam$gene.exposure
  dnam$id.exposure <- dnam$gene.exposure
  
  
  if(psychopathologies[1] !="all")
  {psych <- psych[psych$psych %in% psychopathologies,]}
  
  for(i in 1:nrow(psych)){
    dat <- psych[i,]
    cat("Beginning analysis for:", as.character(dat$psych), "\n")
    if(dat$note %in% c("Z","Zscore","Z-score")){
      cat("ERROR: cannot currently accomodate Z scores. Please refer to prepareMR script. Skipping for now \n")
      next()
    }
    if(dat$psych %in% skipOutcome){
      cat("Skipping this outcome, as per skipOutcome command \n")
      next()
    }
    
    cat("Reading in the outcome data \n")
    #because the files can be very large, this step can be used to reduce the memory footprint of read_outcome_data
    if(largeFiles == T){
      out_temp <- fread(paste0(path_to_stats, dat$filename))
      out_temp <- data.frame(out_temp)
      colnames(out_temp)[colnames(out_temp) == dat$snp_col] <- "rsid"
      rsid <- c(out_temp[,"rsid"])
      out_temp2 <- out_temp[rsid %in% dnam$SNP ,]
      colnames(out_temp2)[colnames(out_temp2) == "rsid"] <- as.character(dat$snp_col)
      write.table(out_temp2, file ="temporary_outcome.txt", sep = "\t", 
                  quote=F, col.names=T, row.names=F)
      rm(out_temp, out_temp2, rsid)
      
      outcome <- read_outcome_data(snps = dnam$SNP,
                                   #filename= paste0(path_to_stats, dat$filename),
                                   filename = "temporary_outcome.txt",
                                   #sep= dat$sep,
                                   sep ="\t",
                                   snp_col = as.character(dat$snp_col),
                                   beta_col = as.character(dat$beta_col),
                                   se_col = as.character(dat$se_col),
                                   effect_allele_col= as.character(dat$effect_allele_col),
                                   other_allele_col = as.character(dat$other_allele_col),
                                   eaf_col = as.character(dat$eaf_col),
                                   pval_col = as.character(dat$pval_col))
    #if(multiVariate==F){ #keep if doing multivariate
      file.remove("temporary_outcome.txt") #clean up the temporary outcome file 
    #}
    }
    else{
      outcome <- read_outcome_data(snps = dnam$SNP,
                                   filename= paste0(path_to_stats, dat$filename),
                                   sep= dat$sep,
                                   snp_col = as.character(dat$snp_col),
                                   beta_col = as.character(dat$beta_col),
                                   se_col = as.character(dat$se_col),
                                   effect_allele_col= as.character(dat$effect_allele_col),
                                   other_allele_col = as.character(dat$other_allele_col),
                                   eaf_col = as.character(dat$eaf_col),
                                   pval_col = as.character(dat$pval_col))
    }
    
    outcome$outcome <- dat$psych
    outcome$id.outcome <- dat$psych
    # test <- fread("YOURPATH/Analysis/Data/Summstats/opi.DEPvEXP_EUR.noAF.tbl")
    # head(test)
    # qnorm(p = 8.94e-7/2)
    # test[test$`P-value` == min(test$`P-value`),]
    
    if(dat$note =="OR"){
      cat("Converting OR to beta \n")
      outcome$beta.outcome <- log(outcome$beta.outcome)
      #outcome$se.outcome <- log(outcome$se.outcome)
    }
    
    cat("Harmonizing data between exposure and outcome \n")
    harm <- harmonise_data(exposure  = dnam, 
                           outcome = outcome)
    if(clumpMR == T){
      harm <- clump_data(harm, clump_r2 = clump.r2)
      #NOTE TO FIX FOR TRANSETHNIC GWAS
    }

    cat("Running MR analyses \n")
    results <- mr(harm)
    resultsSingle <- mr_singlesnp(harm)
    pleiotropy <- mr_pleiotropy_test(harm)
    heterogeneity <- mr_heterogeneity(harm)
    
    if(multiVariate == T){
      cat("Multivariable MR analyses are currently not supported - stay tuned... \n")
    #   cat("Running multivariable MR analyses \n")
    #   #need to get all SNPs for all CpGs 
    #   #
    #   snps <- unique(harm$SNP)
    #   
    #   cgID <- unique(dnam$id.exposure)
    #   cat("Obtaining mQTLs from goDMC - unclumped \n")
    #   query <- list(
    #     rsids = snps,
    #     pval=0.99,
    #     clumped = 0
    #   )
    #   res <- POST("http://api.godmc.org.uk/v0.1/query", body = query, encode = "json")
    #   exposureSNPs <- content(res) %>% lapply(., as_data_frame) %>% bind_rows
    #   dim(exposureSNPs)
    #   
    #   head(exposureSNPs)
    #   mvdat <- mv_harmonise_data(exposure_dat = dnam, outcome_dat = outcome)
    # 
    #   resultsMV <- mv_multiple(mvdat)
    #   
    #   save(resultsMV, 
    #        file = paste0("MVMRresults_", dat$psych, "_", newDir,".Rdata"))
    # 
    }

    cat("Saving results \n \n")
    save(results, resultsSingle, harm, pleiotropy, heterogeneity,
         file = paste0("MRresults_", dat$psych, "_", newDir,".Rdata"))

    rm(outcome, harm, results, resultsSingle, dat)
  }
  
  cat("Done! \n")
  cat(paste0("You can find your results in the following folder: \n",
             paste0(outputdir, newDir), "\n"))
  print(Sys.time()-t)
  sink(type= "output")
  sink(type="message")
  
}
#