### sample script ###
#Loading the dataset
setwd("") #set your working directory!

load("SLCMA_toy_data_2021-03-10.Rdata", verbose=T)
 ## there are two files in this Rdata file, one for the sample data and one for the DNAm data
 #NOTE: ID column must be named "ID" 
head(slcma.sample.data)
### Preparing the covariates ###
colnames(slcma.sample.data)
covars <- c("Sex", "Ethnicity", "Previous.pregnancies", "Maternal.education",
            "Maternal.age", "Birthweight.grams", "Smoked.pregnancy",
            "Bcell", "CD4T", "CD8T", "Gran", "Mono", "NK" #cell types estimated from the Houseman method
)

#any categorical covariate must be structured as a factor
str(slcma.sample.data[,covars])
for(i in c("Sex", "Ethnicity", "Maternal.education","Maternal.age", "Smoked.pregnancy")){
  slcma.sample.data[,i] <- as.factor(slcma.sample.data[,i])
}
str(slcma.sample.data[,covars])


### Preparing the adversities ###
# a list for the recency vectors - there are two adversities in this toy dataset, mompsych and oneadult
colnames(slcma.sample.data)[grep("mompsych", colnames(slcma.sample.data))]
colnames(slcma.sample.data)[grep("oneadult", colnames(slcma.sample.data))]

#make sure this list is in the same order as the grepped column names above
recency.vec.list <- list(c(6, 24, 36, 60)/12, #age in months 
                         c(6, 24, 36, 60)/12)
names(recency.vec.list) <- c("mompsych", "oneadult")



### Preparing the DNA methylation data ###
# outcome vec - vector of CpG site names (~450K)
outcome.vec <- colnames(slcma.DNAm.data[,-1])
# Please NOTE: for this function to work properly, CpG site names must begin with "cg", otherwise, LARS function will need to be adjusted.


### Preparing the data structure for input into the SLCMA ###
#make sure the DNAm data and sample data are in the same order
identical(slcma.sample.data$ID, slcma.DNAm.data$ID) #TRUE

#combine into one data frame - the slcma.df will be used in the SLCMA function
slcma.df <- cbind(slcma.sample.data, slcma.DNAm.data[,-1]) #removing the first column with ID in DNAm data

dim(slcma.df)

### Running the SLCMA ###
adver <- "mompsych" #make sure this is a unique identifier for the columns that are for the exposure (grep is used in the SLCMA)
#adver <- "oneadult" #this would be used in the 

source("LARS-noimpute-function-20210108.R") 
res <- select.LARS.complete(df = slcma.df, 
                            outcome.vec = outcome.vec, 
                            adver = adver,
                            covars = covars,
                            hypos <- c("accumulation", "recency"),
                            exposures <- "default",
                            recency.vec = recency.vec.list[adver][[1]],
                            inf.method = "sI",
                            FWL = TRUE)

#WARNINGS ARE EXPECTED - DON'T WORRY!
#save the file!
save(res, file = paste0("Cohort-SLCMA-", adver, "-",Sys.Date()))

