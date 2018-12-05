
setwd("~/NDN_FASD_I/DMRcate/")
load("~/NDN_FASD_I/Betas.RData")

load("~/NDN_FASD_I/final2.RData")

betas(final2)
dim(pData(final2))
write.table(pData(final2), "~/FASD/filtered.samples.data.txt", sep="\t")

source("http://bioconductor.org/biocLite.R")
biocLite("minfi")
library(DMRcate)
library(Gviz)
library(plyr)
library(methylumi)


identical(betas(final2), (Betas)) #FALSE
identical(rownames(betas(final2)), rownames(Betas)) #TRUE
identical(colnames(betas(final2)), colnames(Betas)) #TRUE

Betas.corr <- betas(final2)
head(Betas.corr)
str(final2)
(featureData(final2))[[5]]
str(featureData(final2))


#transform betas into M-Values
elodie_MS <- logit2(Betas.corr) 
dim(elodie_MS) # 206 samples, 404,030 probes
groups <- read.table("~/NDN_FASD_I/Sample_Groups.txt", header=F)
colnames(groups) <- c("samples", "group")
head(groups)


elodie_MS <- elodie_MS[,order(colnames(elodie_MS))]
colnames(elodie_MS)
groups <- groups[order(as.character(groups$sample)),]
head(groups)
identical(colnames(elodie_MS), as.character(groups$sample))
elodie.design <- model.matrix(~group, data=groups)
load("~/NDN_FASD_I/svs2.RData")
head(svs2)
svs2 <- svs2[order(rownames(svs2)),]
identical(colnames(elodie_MS), rownames(svs2))
elodie.design.2 <- cbind(elodie.design, svs2)
head(elodie.design.2)

elodie.myannotation <- cpg.annotate(datatype = "array", elodie_MS, analysis.type="differential", design=elodie.design.2, coef=2)
elodie.dmrcoutput <- dmrcate(elodie.myannotation, lambda=1000, C=2, pcutoff =0.05)
str(elodie.dmrcoutput)
elodie.results.ranges <- extractRanges(elodie.dmrcoutput, genome = "hg19")
length(elodie.results.ranges) #5702 @ p< 0.05
two.results <- elodie.results.ranges[which(elodie.results.ranges$no.cpgs>1),]
length(two.results) #3005 with two or more probes 


### MORE ETHNICALLY-HOMOGENOUS SUBGROUP ###

cluster1 <- read.table("~/NDN_FASD_I/Cluster1_Samples.txt")
subgroup <- groups[which(groups$sample %in% cluster1$V1),]
head(subgroup)

sub.betas <- Betas.corr[,which(colnames(Betas.corr) %in% subgroup$sample)]
dim(sub.betas)
sub.mval <- logit2(sub.betas)
dim(sub.mval)
sub.mval <- sub.mval[,order(colnames(sub.mval))]
colnames(sub.mval)
subgroup <- subgroup[order(as.character(subgroup$sample)),]
identical(colnames(sub.mval), as.character(subgroup$sample))

sub.elodie.design <- model.matrix(~group, data=subgroup)
load("~/NDN_FASD_I/svs.cauc.RData")
svs.cauc <- svs.cauc[order(rownames(svs.cauc)),]
identical(colnames(sub.mval), rownames(svs.cauc))
sub.elodie.design.2 <- cbind(sub.elodie.design, svs.cauc)
head(sub.elodie.design.2)
library(DMRcate)
sub.elodie.myannotation <- cpg.annotate(datatype = "array",sub.mval, analysis.type="differential", design=sub.elodie.design.2, coef=2)
sub.elodie.dmrcoutput <- dmrcate(sub.elodie.myannotation, lambda=1000, C=2, pcutoff =0.05)
sub.elodie.results.ranges <- extractRanges(sub.elodie.dmrcoutput, genome = "hg19")
length(sub.elodie.results.ranges) #406 @ p<0.05 
two.sub.results <- sub.elodie.results.ranges[which(sub.elodie.results.ranges$no.cpgs>1),]
length(two.sub.results) #289 
unique(two.sub.results$gene_assoc)

length(which(ranges(sub.elodie.results.ranges) %in% ranges(elodie.results.ranges))) #358 @ p<0.05
sig.two <- two.results[(which(ranges(two.results) %in% ranges(two.sub.results))),]
length(sig.two) #101
unique(sig.two$gene_assoc)

sig.two[which(duplicated(sig.two$gene_assoc)),]$gene_assoc
summary(sig.two$no.cpgs) # 2-20 probes, mean = 4.99, median = 4
quantile(sig.two$no.cpgs) # 50% with 4 or more, 25% with 6 or more
summary(sig.two$maxbetafc)
summary(width(ranges(sig.two))) ##41-2450bp, mean = 471bp


RSobject <- RatioSet(Betas.corr, annotation =c(array = "IlluminaHumanMethylation450k", annotation = "ilmn12.hg19"))
RSanno <- getAnnotation(RSobject)

sig.cpgs <- lapply(1:length(sig.two), function(x){
  print(x)
  coords <- names(ranges(sig.two))[x]
  chr <- sub(":.*", "", coords) 
  bookends <- sub(".*:", "", coords)
  startcpg <- as.integer(sub("-.*", "", bookends))
  stopcpg <- as.integer(sub(".*-", "", bookends))
  cpgs <- rownames(RSanno)[RSanno$chr %in% chr & RSanno$pos >= 
                             startcpg & RSanno$pos <= stopcpg]
})

length(unlist(sig.cpgs)) #1748 total CpGs
length(unique(unlist(sig.cpgs))) #980 unique CpGs
sig.cpgs.list <- unlist(sig.cpgs)
DM.probes <- read.table("~/NDN_FASD_I//Significant_Probes.txt")
length(which(sig.cpgs.list %in% DM.probes$V1)) # 91 probes in DMRs that were also in the DM probes

RSanno.sub <- rownames(RS)
sig.annot <- data.frame(cpg = RSanno$Name[which(rownames(RSanno) %in% sig.cpgs.list)],
                        chr= RSanno$chr[which(rownames(RSanno) %in% sig.cpgs.list)],
                        pos =RSanno$pos[which(rownames(RSanno) %in% sig.cpgs.list)])
head(sig.annot)
sig.cpgs.fc <- data.frame(cpg=elodie.dmrcoutput$input$ID[which(elodie.dmrcoutput$input$ID %in% sig.cpgs.list)],
                          FC = elodie.dmrcoutput$input$betafc[which(elodie.dmrcoutput$input$ID %in% sig.cpgs.list)])
head(sig.cpgs.fc)

sig.cpgs.fc <- merge(sig.cpgs.fc, sig.annot, by="cpg")
head(sig.cpgs.fc)
sig.cpgs.fc <- sig.cpgs.fc[order(match(sig.cpgs.fc$cpg, sig.cpgs.list)),]
identical(as.character(sig.cpgs.fc[order(match(sig.cpgs.fc$cpg, sig.cpgs.list)),]$cpg), sig.cpgs.list)


sig.cpg.dmr <- lapply(1:nrow(sig.cpgs.fc), function(x){
  print(x)
  a<- which(sub(":.*", "", names(ranges(sig.two))) == sig.cpgs.fc$chr[x] & 
              start(ranges(sig.two)) <= sig.cpgs.fc$pos[x] & 
              end(ranges(sig.two)) >= sig.cpgs.fc$pos[x])
  names(ranges(sig.two))[a]
})
sig.cpgs.fc$dmr <- unlist(sig.cpg.dmr)
head(sig.cpgs.fc)

dmr.mean.fc <- lapply(1:length(unique(sig.cpgs.fc$dmr)), function(x) {
  mean(sig.cpgs.fc$FC[which(sig.cpgs.fc$dmr == unique(sig.cpgs.fc$dmr)[x])])
})

dmr.fc <- data.frame(dmr = unique(sig.cpgs.fc$dmr),
                     mean.fc = unlist(dmr.mean.fc))

head(dmr.fc)
dmr.fc <- dmr.fc[order(match(dmr.fc$dmr, names(ranges(sig.two)))),]
identical(as.character(dmr.fc$dmr), names(ranges(sig.two)))


sig.two.table <- data.frame(seq = names(ranges(sig.two)),
                            chr = sub(":.*", "", names(ranges(sig.two))),
                            start = start(ranges(sig.two)),
                            end = end(ranges(sig.two)),
                            width = width(ranges(sig.two)),
                            no.probes = sig.two$no.probes,
                            gene_assoc = sig.two$gene_assoc,
                            location = sig.two$group,
                            min_FDR = sig.two$minpval,
                            mean_FDR = sig.two$meanpval,
                            mean_beta_fc = dmr.fc$mean.fc,
                            max_beta_FC = sig.two$maxbetafc)

head(sig.two.table)
dim(sig.two.table)
unique(sig.two.table$location)
write.table(sig.two.table, "~/FASD/DMRcate/significant.DMRs.corrected.txt",sep="\t")

hits <- read.table("~/NDN_FASD_I/DMRcate/significant.DMRs.corrected.txt", sep="\t", header=T)
head(hits)

length(grep("TSS", sig.two.table$location)) #29 TSS 
length(grep("TSS200", sig.two.table$location)) #20 TSS200
length(grep("TSS1500", sig.two.table$location)) #25 TSS1500
length(grep("Body", sig.two.table$location)) #49 body
length(grep("1stExon", sig.two.table$location)) #16 1st exon
length(grep("5'UTR", sig.two.table$location)) #23 5'UTR
length(grep("3'UTR", sig.two.table$location)) #6 3'UTR
length(which(sig.two.table$location=="")) #27 intergenic



source("~/DMRcate_scripts/plot.DMR.R")
plot.DMR(two.sub.results, which(two.sub.results$gene_assoc == "HLA-DPB1"), groups = groups, betas=Betas.corr, toscale=T)
plot.DMR(two.sub.results, grep("NOS1AP", two.sub.results$gene_assoc), groups = groups, betas=Betas.corr, toscale=T)
plot.DMR(two.sub.results, grep("ITGAL", two.sub.results$gene_assoc), groups = groups, betas=Betas.corr, toscale=T) ###
plot.DMR(two.sub.results, grep("UCN3", two.sub.results$gene_assoc), groups = groups, betas=Betas.corr, toscale=T) ###

plot.DMR(sig.two, which(sig.two$gene_assoc == "SLC38A2"), groups = groups, betas=Betas.corr, toscale=T, many.genes=F)
plot.DMR(sig.two, which(sig.two$gene_assoc == "RNMTL1"), groups = groups, betas=Betas.corr, toscale=T)
plot.DMR(sig.two, which(sig.two$gene_assoc == "SNED1"), groups = groups, betas=Betas.corr, toscale=T)
plot.DMR(sig.two, grep("PCDH12", sig.two$gene_assoc), groups = groups, betas=Betas.corr, toscale=T)
plot.DMR(sig.two, grep("PCDHA", sig.two$gene_assoc), groups = groups, betas=Betas.corr, toscale=T)
plot.DMR(sig.two, grep("SLC22A18", sig.two$gene_assoc), groups = groups, betas=Betas.corr, toscale=T)
plot.DMR(sig.two, grep("HKR1", sig.two$gene_assoc), groups = groups, betas=Betas.corr, toscale=T)
plot.DMR(sig.two, grep("WDR52", sig.two$gene_assoc), groups = groups, betas=Betas.corr, toscale=T)
plot.DMR(sig.two, grep("KCNAB2", sig.two$gene_assoc), groups = groups, betas=Betas.corr, toscale=T)
plot.DMR(sig.two, grep("KCNQ1", sig.two$gene_assoc), groups = groups, betas=Betas.corr, toscale=T)
plot.DMR(sig.two, grep("NECAB3", sig.two$gene_assoc), groups = groups, betas=Betas.corr, toscale=T)
plot.DMR(sig.two, grep("PRKDC", sig.two$gene_assoc), groups = groups, betas=Betas.corr, toscale=T)
plot.DMR(sig.two, grep("STRA6", sig.two$gene_assoc), groups = groups, betas=Betas.corr, toscale=T) 
plot.DMR(sig.two, grep("HEATR2", sig.two$gene_assoc), groups = groups, betas=Betas.corr, toscale=T)
plot.DMR(sig.two, grep("IFT140", sig.two$gene_assoc), groups = groups, betas=Betas.corr, toscale=T)
plot.DMR(sig.two, grep("RADIL", sig.two$gene_assoc), groups = groups, betas=Betas.corr, toscale=T)
plot.DMR(sig.two, grep("HKDC1", sig.two$gene_assoc), groups = groups, betas=Betas.corr, toscale=T)
plot.DMR(sig.two, grep("ARHGEF19", sig.two$gene_assoc), groups = groups, betas=Betas.corr, toscale=T)
plot.DMR(sig.two, grep("NDST4", sig.two$gene_assoc), groups = groups, betas=Betas.corr, toscale=T)
plot.DMR(sig.two, grep("EID3", sig.two$gene_assoc), groups = groups, betas=Betas.corr, toscale=T)
plot.DMR(sig.two, grep("RGL3", sig.two$gene_assoc), groups = groups, betas=Betas.corr, toscale=T)
plot.DMR(sig.two, grep("PPP1R2P1", sig.two$gene_assoc), groups = groups, betas=Betas.corr, toscale=T)


sig.two.table[which(duplicated(sig.two.table$gene_assoc)),]$gene_assoc
source("~/FASD/DMRcate/plot.chr.R")
sig.two.table[grep("IFT140,TMEM204", sig.two$gene_assoc),] #12 probes total
plot.DMR(sig.two, grep("IFT140,TMEM204", sig.two$gene_assoc)[1], groups = groups, betas=Betas.corr, toscale=T)
plot.DMR(sig.two, grep("IFT140,TMEM204", sig.two$gene_assoc)[2], groups = groups, betas=Betas.corr, toscale=T)
plot.chr(chrom=16, start=1583810, end = 1599150, groups = groups, betas=Betas.corr, toscale=T)
sig.two.table[grep("UPP1", sig.two$gene_assoc),] #6 probes total
plot.chr(chrom=7, start=48129814, end = 48132109, groups = groups, betas=Betas.corr)
sig.two.table[grep("RNMTL1", sig.two$gene_assoc),] #7 probes total
plot.chr(chrom=17, start=693137, end = 695661, groups = groups, betas=Betas.corr, toscale=T) #REALLY COOL! 2 DMRs with opposite directions
