### Making the sex-diff peakset ###
library(DiffBind)
cat("Preparing index for sex-differences analysis", '\n')

index <- read.table("location_of_data/dev2_seq/index.corr.txt", sep='\t', header=T)
diff.index <- index[which(index$tissue %in% c("PFC")),]
dim(diff.index) #30,11
diff.index <- diff.index[order(diff.index$sex, diff.index$group, diff.index$sampID),]
diff.index$replicate <- c(rep(1:5,6))

for(i in 1:nrow(diff.index)){
  cat(paste("Now working on sample #", i, ": ", diff.index$sampID[i], " ", diff.index$sex[i], sep="" ), "\n")
  
  if(i ==1){diff.peakset <- dba.peakset(peaks = as.character(diff.index$peaks[i]), sampID = as.character(diff.index$sampID[i]), 
                                        tissue= as.character(diff.index$tissue[i]), treatment= as.character(diff.index$group[i]), 
                                        replicate= as.numeric(diff.index$replicate[i]), condition = as.character(diff.index$sex[i]),
                                        peak.caller="macs", peak.format="macs", bamReads= as.character(diff.index$reads[i]))
  }
  
  diff.peakset <- dba.peakset(diff.peakset, peaks = as.character(diff.index$peaks[i]), sampID = as.character(diff.index$sampID[i]), 
                              tissue= as.character(diff.index$tissue[i]), treatment= as.character(diff.index$group[i]), 
                              replicate= as.numeric(diff.index$replicate[i]), condition = as.character(diff.index$sex[i]),
                              peak.caller="macs", peak.format="macs", bamReads= as.character(diff.index$reads[i]))
}

cat("Sex differences peakset successfully generated!","\n")
cat("Saving sex differences peakset","\n")

save(diff.peakset, file="location_of_data/dev2_seq/dev2.diff.peakset.Rdata")

cat("Selecting peaks with minimum overlap of 3 samples","\n")
diff.peakset.OL3 = dba.count(diff.peakset, minOverlap = 3)

cat("Saving final male peakset and index")
save(diff.peakset.OL3, file="location_of_data/dev2_seq/dev2.diff.peakset.OL3.Rdata")
write.table(diff.index, file = "location_of_data/dev2_seq/dev2.diff.index.txt", sep='\t', quote=F, col.names = T)

cat("Cleaning up","\n")
rm(diff.peakset)
rm(diff.peakset.OL3)
rm(diff.index)

cat("Done!","\n")

