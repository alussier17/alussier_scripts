
#compile bias results

setwd("/data/js95/DunnLab/alussier/DCHS_technical/bias_analysis/results")

files <- list.files(pattern = ".Rdata")

for(i in files){
  if(i == files[1]){biasResults <- data.frame(); biasTukey <- data.frame(); failed <- c()}
  load(i)
  if(nrow(anovaSummary)==0){
    failed <- c(failed, i)
  }
  biasResults <- rbind(biasResults, anovaSummary)
  biasTukey <- rbind(biasTukey, anovaTukey)
}
dim(biasResults)
dim(biasTukey)
length(failed)
failed


#checking the failed CpGs - could be that they are unique to each array
colSums(is.na(fullBeta2))

load("../betasFull_2023-10-18.Rdata")
table(rowSums(is.na(fullBeta)))[3] #number with only one array
nrow(fullBeta) - 252523 - nrow(biasResults) #all good

save(biasResults, biasTukey, 
     file="/data/js95/DunnLab/alussier/DCHS_technical/bias_analysis/biasAnalysisResults_2023-10-26.Rdata")
list.files()

