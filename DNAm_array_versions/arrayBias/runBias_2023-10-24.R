
runBias <- function(start, end){
  cat("Starting setup at: ", paste0(Sys.time()), "\n")
  library(broom)
  library(dplyr)
  library(furrr)
  library(limma)
  library(purrr)
  library(reshape)
  library(tidyr)
  
  setwd("/data/js95/DunnLab/alussier/DCHS_technical/bias_analysis/")
  load("../betasFull_2023-10-18.Rdata")
  
  allCpGs <- fullBeta$rowname
  if(end > nrow(fullBeta)){end <- nrow(fullBeta)}
  allCpGs2 <- allCpGs[start:end]  
  
  fullBeta2 <- fullBeta[start:end,]
  rm(fullBeta)
  
  t <- Sys.time()
  cat("Starting main analysis at: ", paste0(Sys.time()), "\n")
  results <- future_map(allCpGs2, ~{
    x = .x
    dat <- fullBeta2[fullBeta2$rowname == x, ]
    dat <- suppressMessages(melt(dat))
    dat$Array <- as.factor(strsplit2(dat$variable, "_")[, 2])
    dat$SampleID <- as.factor(strsplit2(dat$variable, "_")[, 1])
    dat <- na.omit(dat)
    if (length(unique(dat$Array)) == 1) {
      return(NULL)
    }
    a <- aov(value ~ Array + SampleID, data = dat)
    b <- broom::tidy(a)
    b$CpG <- x
    
    c <- data.frame(TukeyHSD(a)$Array)
    c$CpG <- x 
    c$Contrast <- rownames(c)
    
    list(anova = b[1, ], tukey = c)
  })
  print(Sys.time()-t)
  
  cat("Combining results at: ", paste0(Sys.time()), "\n")
  # Combine the results
  anovaSummary <- data.frame()
  anovaTukey <- data.frame()
  
  for (i in seq_along(allCpGs2)) {
    if (!is.null(results[[i]]$anova)) {
      anovaSummary <- rbind(anovaSummary, results[[i]]$anova)
    }
    if (!is.null(results[[i]]$tukey)) {
      anovaTukey <- rbind(anovaTukey, results[[i]]$tukey)
    }
  }
  
  dim(anovaSummary)
  dim(anovaTukey)
  
  save(anovaSummary, anovaTukey, file= paste0("results/anovaBias_", start,"_",end, ".Rdata"))
  cat("Done at: ", paste0(Sys.time()), "\n")
}

