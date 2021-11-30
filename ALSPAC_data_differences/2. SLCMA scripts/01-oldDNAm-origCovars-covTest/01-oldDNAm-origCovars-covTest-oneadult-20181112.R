
setwd("/data/js95/ALSPAC/ARIES/F7/BPrev2")
## packages 

## data
load("ARIES_betas_noTwin_rfaminst_pcs_20181023.Rdata")
outcome.vec <- colnames(df)[87:ncol(df)]

## load function for stage 1 LARS, complete data
source("LARS-noimpute-function-20181023.R")

adver <- "oneadult"

## a list for the recency vectors 
recency.vec.list <- list(c(18, 30, 42, 57, 69, 81)/12,
                   c(8, 21, 33, 61, 85)/12,
                   c(8, 21, 33, 47, 85)/12,
                   c(18, 30, 42, 57, 69, 81)/12,
                   c(21, 33, 61, 85)/12,
                   c(8, 21, 33, 47, 61, 73)/12,
                   c(8, 21, 33, 61, 73)/12)
names(recency.vec.list) <-   c("abuse", "Fscore",  "oneadult", "r_faminst", "nbhqual","parc","mompsy")


res <- select.LARS.complete(df, 
                                  outcome.vec = outcome.vec, 
                                  adver = adver,
                                  covars <- c(## baseline characteristics,EXCLUDE SES
                                    "SES_parent", "WHITE", "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
                                    # cell counts
                                    "CD8", "CD4", "CD56", "CD19", "CD14", "Gran"),
                                  hypos <- c("accumulation", "recency"),
                                  exposures <- "default",
                                  recency.vec = recency.vec.list[adver][[1]],
                                  inf.method = "covTest")

save(res, file = paste0("LARS_rev2_rerun_stage1_covTest_", adver, "_20181112.Rdata"))

