### MAKING A TOY DATASET FOR SLCMA ###
setwd("~/Dropbox (Partners HealthCare)/Dunn_shared/PACE Consortium/SLCMA_sample")
set.seed(17)
slcma.sample.data <- data.frame(ID =  paste0(sample(1:10000,100), sample(rep(c("A","A","A","B"),25), 100)))
head(slcma.sample.data)

#Fake covariates
slcma.sample.data$Sex <- sample(rep(c(0,1),50), 100)
slcma.sample.data$Ethnicity <- sample(rep(c(1,1,1,0),25 ), 100)
slcma.sample.data$Previous.pregnancies <- sample(rep(c(1,1,2,3,4), 20), 100)
slcma.sample.data$Maternal.education <- sample(rep(c(0,1,2,3),25), 100)
slcma.sample.data$Maternal.age <- c(sample(rep(c(0,1,2),33), 99),1)
slcma.sample.data$Birthweight.grams <- round(runif(100, 2500, 4000),0)
slcma.sample.data$Smoked.pregnancy <- sample(rep(c(1,1,1,0),25 ), 100)

#random cell types, taken from F7, but scrambled to at least have real distributions
cells.f7 <- read.table("~/Dropbox (Partners HealthCare)/Dunn_shared/ALSPAC_DNAm_differences/data/celltypes_F7_20200121.txt",
                       header=T)
slcma.sample.data <- cbind(slcma.sample.data, 
                           cells.f7[sample(1:nrow(cells.f7),100), -1])
rm(cells.f7)

head(slcma.sample.data)

#fake outcome data for two adversities
slcma.sample.data$mompsych.6m <- sample(rep(c(0,0,0,0,1),20), 100)
slcma.sample.data$mompsych.24m <- sample(rep(c(0,0,0,0,1),20), 100)
slcma.sample.data$mompsych.36m <- sample(rep(c(0,0,0,0,1),20), 100)
slcma.sample.data$mompsych.60m <- sample(rep(c(0,0,0,0,1),20), 100)

slcma.sample.data$oneadult.6m <- sample(rep(c(0,0,0,0,1),20), 100)
slcma.sample.data$oneadult.24m <- sample(rep(c(0,0,0,0,1),20), 100)
slcma.sample.data$oneadult.36m <- sample(rep(c(0,0,0,0,1),20), 100)
slcma.sample.data$oneadult.60m <- sample(rep(c(0,0,0,0,1),20), 100)

dim(slcma.sample.data)
head(slcma.sample.data)

#making fake CpG data
for(i in 1:100){
  if(i==1){slcma.DNAm.data <- data.frame(ID= slcma.sample.data$ID)}
  set.seed(i)
  min = sample(seq(0,0.95, by= 0.05),1)
  max = sample(seq(0.05,1, by= 0.05),1)
  if(min == max){max <- min + 0.05}
  slcma.DNAm.data$temp <- runif(100, min = ifelse(min<max, min, max), max = ifelse(min<max, max, min))
  colnames(slcma.DNAm.data)[1+i] <- ifelse(i<10, paste0("cg00",i), 
                                           ifelse(i<100, paste0("cg0",i), paste0("cg", i)))
}
summary(slcma.DNAm.data[,-1])


#saving as a single dataset
save(slcma.sample.data, slcma.DNAm.data, file= "SLCMA_toy_data_2021-03-10.Rdata")
