
setwd("/data/js95/ALSPAC/ARIES/F7/biopsy2019-rerun/")
## packages 

## data
# load("ARIES_betas_noTwin_rfaminst_pcs_20181023.Rdata")
load("ARIES_BP_rerun_20200225.Rdata")

# colnames(aries)

outcome.vec <- colnames(aries)[102:ncol(aries)]

## double check the variable types
# lapply(aries[,1:101], function(x){
#   class(x)
# })

## load function for stage 1 LARS, complete data
source("LARS-noimpute-function-20200225.R")

adver <- "abuse"

## a list for the recency vectors 
recency.vec.list <- list(c(18, 30, 42, 57, 69, 81)/12,
                   c(8, 21, 33, 61, 85)/12,
                   c(8, 21, 33, 47, 85)/12,
                   c(18, 30, 42, 57, 69, 81)/12,
                   c(21, 33, 61, 85)/12,
                   c(8, 21, 33, 47, 61, 73)/12,
                   c(8, 21, 33, 61, 73)/12)
names(recency.vec.list) <-   c("abuse", "Fscore",  "oneadult", "r_faminst", "nbhqual","parc","mompsy")


res <- select.LARS.complete(df = aries, 
                                  outcome.vec = outcome.vec, 
                                  adver = adver,
                                  covars = c(## baseline characteristics,EXCLUDE SES
                                    "SES_parent", "WHITE", "Female", "mom_birthage", "ppregnum", "birth.weight","sustained.smoke",
                                    # cell counts
                                    "Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK"),
                                  hypos <- c("accumulation", "recency"),
                                  exposures <- "default",
                                  recency.vec = recency.vec.list[adver][[1]],
                                  inf.method = "covTest")

save(res, file = paste0("SLCMA-newDNAm-rerun-", adver, "_20200225.Rdata"))

