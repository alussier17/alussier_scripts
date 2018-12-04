#MAKING SOME ANNOTATION - P22 #
seqs.P22 <- as.matrix(paste(P22H_WBC.peakset$peaks[[1]]$Chr, P22H_WBC.peakset$peaks[[1]]$Start, 
                            P22H_WBC.peakset$peaks[[1]]$End, sep=":"))
seqs.P22 <- cbind(as.character(P22H_WBC.peakset$peaks[[1]]$Chr), as.numeric(P22H_WBC.peakset$peaks[[1]]$Start), 
                  as.numeric(P22H_WBC.peakset$peaks[[1]]$End), seqs.P22)
colnames(seqs.P22) <- c("chr", "start","end", "genes")
head(seqs.P22)
write.table(seqs.P22, file = "~/region_bams/P22_peaks_annot/P22.peakset.bed", sep="\t", row.names=F, quote=F, col.names=F)
seqs.P22 <- read.table("~/region_bams/P22_peaks_annot/P22.peakset.bed")
colnames(seqs.P22) <- c("chr", "start","end", "genes")

P22.introns <- read.table("~/region_bams/P22_peaks_annot/P22.peakset.introns.bed",sep="\t", header=F)
P22.introns <- P22.introns[,c(-1,-2,-3,-5,-6,-7,-9,-10)]
colnames(P22.introns) <- c("genes","introns")
dim(P22.introns)
head(unique(P22.introns))
introns.agg <- aggregate(introns~genes, data=P22.introns, FUN=paste, collapse=" | ")
head(introns.agg)
dim(introns.agg)

P22.exons  <- read.table("~/region_bams/P22_peaks_annot/P22.peakset.exons.bed",sep="\t", header=F)
P22.exons <- P22.exons[,c(-1,-2,-3,-5,-6,-7,-9,-10)]
colnames(P22.exons) <- c("genes","exons")
dim(P22.exons)
head(P22.exons)
exons.agg <- aggregate(exons~genes, data=P22.exons, FUN=paste, collapse=" | ")
head(exons.agg)

P22.cpg <- read.table("~/region_bams/P22_peaks_annot/P22.peakset.cpg.bed",sep="\t", header=F)
P22.cpg <- cbind(P22.cpg, paste(P22.cpg[,5], P22.cpg[,6], P22.cpg[,7], P22.cpg[,8], sep=":"))
P22.cpg <- P22.cpg[,c(-1,-2,-3,-5,-6,-7,-8)]
colnames(P22.cpg) <- c("genes","cpg")
dim(P22.cpg)
cpg.agg <- aggregate(cpg~genes, data=P22.cpg, FUN=paste, collapse=" | ")
head(cpg.agg)

P22.promoter  <- read.table("~/region_bams/P22_peaks_annot/P22.peakset.prom200.bed",sep="\t", header=F)
head(P22.promoter)
P22.promoter <- P22.promoter[,c(-1,-2,-3,-5,-6,-7,-9,-10)]
colnames(P22.promoter) <- c("genes", "prom200")
P22.promoter$prom200 <- unlist(lapply(1:nrow(P22.promoter), function(x){
  a <- strsplit(as.character(P22.promoter$prom200[x]), "_")[[1]][1:2]
  paste(a[1], a[2], sep="_")}))
head(P22.promoter)
dim(P22.promoter)
prom.agg <- aggregate(prom200~genes, data=P22.promoter, FUN=paste, collapse=" | ")
head(prom.agg)


P22.prom1500  <- read.table("~/region_bams/P22_peaks_annot/P22.peakset.prom1500.bed",sep="\t", header=F)
head(P22.prom1500)
P22.prom1500 <- P22.prom1500[,c(-1,-2,-3,-5,-6,-7,-9,-10)]
colnames(P22.prom1500) <- c("genes", "prom1500")
P22.prom1500$prom1500 <- unlist(lapply(1:nrow(P22.prom1500), function(x){
  a <- strsplit(as.character(P22.prom1500$prom1500[x]), "_")[[1]][1:2]
  paste(a[1], a[2], sep="_")
}))
head(P22.prom1500)
dim(P22.prom1500)
prom1500.agg <- aggregate(prom1500~genes, data=P22.prom1500, FUN=paste, collapse=" | ")
head(prom1500.agg)
dim(prom1500.agg)

P22.three <- read.table("~/region_bams/P22_peaks_annot/P22.peakset.UTR.3.bed",sep="\t", header=F)
head(P22.three)
P22.three <- P22.three[,c(-1,-2,-3,-5,-6,-7,-9,-10)]
colnames(P22.three) <- c("genes","UTR.3")
dim(P22.three)
head(P22.three)
three.agg <- aggregate(UTR.3~genes, data=P22.three, FUN=paste, collapse=" | ")
head(three.agg)
dim(three.agg)

P22.five <- read.table("~/region_bams/P22_peaks_annot/P22.peakset.UTR.5.bed",sep="\t", header=F)
head(P22.five)
P22.five <- P22.five[,c(-1,-2,-3,-5,-6,-7,-9,-10)]
colnames(P22.five) <- c("genes","UTR.5")
dim(P22.five)
five.agg <- aggregate(UTR.5~genes, data=P22.five, FUN=paste, collapse=" | ")
head(five.agg)
dim(five.agg)

head(P22.genes)
P22.genes <- read.table("~/region_bams/P22_peaks_annot/P22.peakset.genes.bed",sep="\t", header=F)
P22.genes <- P22.genes[,c(-1,-2,-3,-5,-6,-7,-9,-10,-11,-12,-13,-14,-15,-16)]
colnames(P22.genes) <- c("genes","ID")
colnames(P22.promoter) = colnames(P22.genes)
P22.genes <- rbind(P22.genes, P22.promoter)
dim(P22.genes)
head(P22.genes)


P22.refseq <- read.table("~/region_bams/rn6.refseq.bed",sep="\t", header=F)
P22.refseq <- P22.refseq[,c(-2,-3)]
colnames(P22.refseq) <- c("ID","name")
dim(P22.refseq)
head(P22.refseq)

P22.genes.ID <- merge(P22.genes, P22.refseq, by="ID") 
P22.genes.ID <- P22.genes.ID[-which(duplicated(P22.genes.ID)),]
head(P22.genes.ID)
dim(P22.genes.ID)
P22.genes.ID.agg <- aggregate(ID~genes, data = P22.genes.ID, FUN = paste, collapse = " | ")



P22.regions <- data.frame(genes=seqs.P22[,4])
dim(P22.regions)
head(P22.regions)

P22.annot <- merge(P22.regions, P22.genes.ID.agg, by="genes",all.x=T)
dim(P22.annot)
head(P22.annot)

P22.annot <- merge(P22.annot, prom.agg, by="genes",all=T)
dim(P22.annot)
P22.annot <- merge(P22.annot, prom1500.agg, by="genes",all=T)
dim(P22.annot)
P22.annot <- merge(P22.annot, cpg.agg, by="genes", all=T)
dim(P22.annot)
P22.annot <- merge(P22.annot, exons.agg, by="genes",all=T)
dim(P22.annot)
P22.annot <- merge(P22.annot, introns.agg, by="genes",all.x=T)
dim(P22.annot)
P22.annot <- merge(P22.annot, five.agg, by="genes",all.x=T)
dim(P22.annot)
P22.annot <- merge(P22.annot, three.agg, by="genes",all.x=T)
dim(P22.annot)
P22.annot$exon.1 <- NA
P22.annot$exon.1[grep("exon_1_", P22.annot$exons)] <- "yes"
dim(P22.annot)
head(P22.annot)
head(P22.annot[which(!is.na(P22.annot$exon.1)),])
length(unique(P22.annot$genes))
save(P22.annot, file= "~/region_bams/P22_peaks_annot/P22.peakset.annotation.Rdata")


# P22 annotation of repeats - LINE/SINE/LTR/ETC.... 
repeats <- read.table("~/region_bams/repeat_masker.txt")
head(repeats)
repeats.simple <- repeats[,c(6,7,8,11,12,13)]
head(repeats.simple)
write.table(repeats.simple, file="~/region_bams/repeats.bed", col.names = F, quote=F, row.names=F, sep="\t")
#intersectBed -wa -wb -a ~/region_bams/P22_peaks_annot/P22.peakset.bed -b ~/region_bams/repeats.bed > ~/region_bams/P22_peaks_annot/rn6.repeats.bed
P22.repeats <- read.table("~/region_bams/P22_peaks_annot/rn6.repeats.bed")
head(P22.repeats)
colnames(P22.repeats) <- c("chr","start","stop","genes","chr_rep","start_rep","stop_rep","repName","repClass","repFamily")
P22.repeats <- P22.repeats[,c(4,8:10)]

p22.repeats.agg.name <- aggregate(repName~genes, data=P22.repeats, FUN=paste, collapse=" | ")
p22.repeats.agg.class <- aggregate(repClass~genes, data=P22.repeats, FUN=paste, collapse=" | ")
p22.repeats.agg.family <- aggregate(repFamily~genes, data=P22.repeats, FUN=paste, collapse=" | ")
dim(p22.repeats.agg.name)
dim(p22.repeats.agg.class)
dim(p22.repeats.agg.family)

p22.repeats.agg <- cbind(p22.repeats.agg.name, p22.repeats.agg.class$repClass, p22.repeats.agg.family$repFamily)
colnames(p22.repeats.agg)[3:4] <- c("repClass","repFamily")
head(p22.repeats.agg)
dim(p22.repeats.agg)

save(p22.repeats.agg, file= "~/region_bams/P22_peaks_annot/P22.peakset.annotation.repeats.Rdata")




## MAKING AN ANNOTATION FOR CG CONTENT OF REGIONS ##
p22.chrs<- (unique(seqs[,1]))
p22.chrs<- gsub("KL","kl",p22.chrs)
p22.chrs<- gsub("chrUn","chrun",p22.chrs)
p22.chrs<- gsub("chrX","chrx",p22.chrs)
p22.chrs<- gsub("chrY","chry",p22.chrs)
p22.chrs<- gsub("AABR","aabr",p22.chrs)

rn6 <- list.files("/home/alussier/Rnor6/fasta", full=T)
rn6.p22.chr<- lapply(1:length(p22.chrs), function(x) rn6[grep(paste(p22.chrs[x], ".fasta",sep=""), rn6)])

p22.rn6.fasta <- lapply(rn6.p22.chr, read.fasta)
str(p22.rn6.fasta[[3]])

seqs.p22.cor <- gsub("KL","kl",seqs)
seqs.p22.cor<- gsub("chrUn","chrun",seqs.p22.cor)
seqs.p22.cor<- gsub("chrX","chrx",seqs.p22.cor)
seqs.p22.cor<- gsub("chrY","chry",seqs.p22.cor)
seqs.p22.cor<- gsub("AABR","aabr",seqs.p22.cor)

p22.chrs.cor <- as.factor(seqs.p22.cor[,1])
p22.chrs.cor

p22.sequence <- sapply(1:nrow(seqs), function(x) p22.rn6.fasta[[p22.chrs.cor[x]]][[1]][seqs[x,2]:seqs[x,3]]) # pull reference sequence for feature
p22.CpG <- sapply(1:nrow(seqs), function(x) str_count(paste(p22.sequence[[x]],collapse=''), "cg")) #total number of CpG in reference
p22.CpG.per.bp <- sapply(1:nrow(seqs), function(x) p22.CpG[x]/length(p22.sequence[[x]])) #CpG / base pair in feature
p22.GC.content <- sapply(1:nrow(seqs), function(x) sum(str_count(paste(p22.sequence[[x]],collapse=''), c("c","g")))) # total number of cytosine or guanine nucleotides
p22.GC.feature <- sapply(1:nrow(seqs), function(x) p22.GC.content[x]/length(p22.sequence[[x]])) # GC content of feature (per base pair)

p22.CG.annot <- data.frame(seqs) 
p22.CG.annot$length <- sapply(1:nrow(seqs), function(x) length(p22.sequence[[x]]))
p22.CG.annot$CpG.number <- p22.CpG
p22.CG.annot$CpG.content <- p22.CpG.per.bp
p22.CG.annot$GC.number <- p22.GC.content
p22.CG.annot$GC.content <- p22.GC.feature

p22.CG.annot <- p22.CG.annot[,c(4,1,2,3,5,6,7,8,9)]
head(p22.CG.annot)
save(p22.CG.annot, file= "~/region_bams/P22_peaks_annot/P22.peakset.CG.annotation.Rdata")



## MAKING A RAW SEQUENCES ANNOTATION ##
p22.sequences <- data.frame(seqs)
p22.sequences <- p22.sequences[,c(4,1,2,3)]
p22.sequences$seq <- sapply(1:nrow(seqs), function(x) paste(p22.sequence[[x]],collapse=''))
head(p22.sequences)

save(p22.sequences, file= "~/region_bams/P22_peaks_annot/P22.peakset.sequences.Rdata")

####################################################################################################################################

#MAKING DEV ANNOTATION - BED FILES OF INTRONS/EXONS/PROMOTERS/CPG/GENES ARE MADE USING "intersectBed -wa -wb -a P22.peakset.bed -b tableFromUCSC.bed"
seqs.dev <- as.matrix(paste(all.peakset$peaks[[1]]$Chr, all.peakset$peaks[[1]]$Start, 
                            all.peakset$peaks[[1]]$End, sep=":"))
seqs.dev <- cbind(as.character(all.peakset$peaks[[1]]$Chr), as.numeric(all.peakset$peaks[[1]]$Start), 
                  as.numeric(all.peakset$peaks[[1]]$End), seqs.dev)
colnames(seqs.dev) <- c("chr", "start","end", "genes")
head(seqs.dev)
write.table(seqs.dev, file = "~/region_bams/dev_peaks_annot/dev.peakset.bed", sep="\t", row.names=F, quote=F, col.names=F)
seqs.dev <- read.table("~/region_bams/dev_peaks_annot/dev.peakset.bed")
colnames(seqs.dev) <- c("chr", "start","end", "genes")

dev.introns <- read.table("~/region_bams/dev_peaks_annot/dev.peakset.introns.bed",sep="\t", header=F)
dev.introns <- dev.introns[,c(-1,-2,-3,-5,-6,-7,-9,-10)]
colnames(dev.introns) <- c("genes","introns")
dim(dev.introns)
head(unique(dev.introns))
introns.agg <- aggregate(introns~genes, data=dev.introns, FUN=paste, collapse=" | ")
head(introns.agg)
dim(introns.agg)

dev.exons  <- read.table("~/region_bams/dev_peaks_annot/dev.peakset.exons.bed",sep="\t", header=F)
dev.exons <- dev.exons[,c(-1,-2,-3,-5,-6,-7,-9,-10)]
colnames(dev.exons) <- c("genes","exons")
dim(dev.exons)
head(dev.exons)
exons.agg <- aggregate(exons~genes, data=dev.exons, FUN=paste, collapse=" | ")
head(exons.agg)

dev.cpg <- read.table("~/region_bams/dev_peaks_annot/dev.peakset.cpg.bed",sep="\t", header=F)
dev.cpg <- cbind(dev.cpg, paste(dev.cpg[,5], dev.cpg[,6], dev.cpg[,7], dev.cpg[,8], sep=":"))
dev.cpg <- dev.cpg[,c(-1,-2,-3,-5,-6,-7,-8)]
colnames(dev.cpg) <- c("genes","cpg")
dim(dev.cpg)
cpg.agg <- aggregate(cpg~genes, data=dev.cpg, FUN=paste, collapse=" | ")
head(cpg.agg)

dev.promoter  <- read.table("~/region_bams/dev_peaks_annot/dev.peakset.prom200.bed",sep="\t", header=F)
dev.promoter <- dev.promoter[,c(-1,-2,-3,-5,-6,-7,-9,-10)]
colnames(dev.promoter) <- c("genes", "prom200")
dev.promoter$prom200 <- unlist(lapply(1:nrow(dev.promoter), function(x){
  a <- strsplit(as.character(dev.promoter$prom200[x]), "_")[[1]][1:2]
  paste(a[1], a[2], sep="_")}))
head(dev.promoter)
dim(dev.promoter)
prom.agg <- aggregate(prom200~genes, data=dev.promoter, FUN=paste, collapse=" | ")
head(prom.agg)
dim(prom.agg)

dev.prom1500  <- read.table("~/region_bams/dev_peaks_annot/dev.peakset.prom1500.bed",sep="\t", header=F)
head(dev.prom1500)
dev.prom1500 <- dev.prom1500[,c(-1,-2,-3,-5,-6,-7,-9,-10)]
colnames(dev.prom1500) <- c("genes", "prom1500")
dev.prom1500$prom1500 <- unlist(lapply(1:nrow(dev.prom1500), function(x){
  a <- strsplit(as.character(dev.prom1500$prom1500[x]), "_")[[1]][1:2]
  paste(a[1], a[2], sep="_")
}))
head(dev.prom1500)
dim(dev.prom1500)
prom1500.agg <- aggregate(prom1500~genes, data=dev.prom1500, FUN=paste, collapse=" | ")
head(prom1500.agg)
dim(prom1500.agg)

dev.three <- read.table("~/region_bams/dev_peaks_annot/dev.peakset.UTR.3.bed",sep="\t", header=F)
dev.three <- dev.three[,c(-1,-2,-3,-5,-6,-7,-9,-10)]
colnames(dev.three) <- c("genes","UTR.3")
dim(dev.three)
head(dev.three)
three.agg <- aggregate(UTR.3~genes, data=dev.three, FUN=paste, collapse=" | ")
head(three.agg)
dim(three.agg)

dev.five <- read.table("~/region_bams/dev_peaks_annot/dev.peakset.UTR.5.bed",sep="\t", header=F)
dev.five <- dev.five[,c(-1,-2,-3,-5,-6,-7,-9,-10)]
colnames(dev.five) <- c("genes","UTR.5")
dim(dev.five)
head(dev.five)
five.agg <- aggregate(UTR.5~genes, data=dev.five, FUN=paste, collapse=" | ")
head(five.agg)
dim(five.agg)


dev.genes <- read.table("~/region_bams/dev_peaks_annot/dev.peakset.genes.bed",sep="\t", header=F)
dev.genes <- dev.genes[,c(-1,-2,-3,-5,-6,-7,-9,-10,-11,-12,-13,-14,-15,-16)]
colnames(dev.genes) <- c("genes","ID")
colnames(dev.promoter) = colnames(dev.genes)
dev.genes <- rbind(dev.genes, dev.promoter)
dev.genes <- dev.genes[-which(duplicated(dev.genes)),]
dim(dev.genes)
head(dev.genes)


dev.refseq <- read.table("~/region_bams/rn6.refseq.bed",sep="\t", header=F)
dev.refseq <- dev.refseq[,c(-2,-3)]
colnames(dev.refseq) <- c("ID","name")
dev.refseq <- dev.refseq[-which(duplicated(dev.refseq)),]
head(dev.refseq)

dev.genes.ID <- merge(dev.genes, dev.refseq, by="ID") 
head(dev.genes.ID)
dim(dev.genes.ID)
dev.genes.ID.agg <- aggregate(ID~genes, data=dev.genes.ID, FUN=paste, collapse=" | ")
dim(dev.genes.ID.agg)


dev.regions <- data.frame(genes=seqs.dev[,4])
dim(dev.regions) #469,498
head(dev.regions)

dev.annot <- merge(dev.regions, dev.genes.ID.agg, by="genes", all.x =T)
dim(dev.annot)
head(dev.annot)
dev.annot <- merge(dev.annot, prom.agg, by="genes",all=T)
dim(dev.annot)
dev.annot <- merge(dev.annot, prom1500.agg, by="genes",all=T)
dim(dev.annot)
dev.annot <- merge(dev.annot, cpg.agg, by="genes", all=T)
dim(dev.annot)
dev.annot <- merge(dev.annot, exons.agg, by="genes",all=T)
dim(dev.annot)
dev.annot <- merge(dev.annot, introns.agg, by="genes",all.x=T)
dim(dev.annot)
dev.annot <- merge(dev.annot, five.agg, by="genes",all.x=T)
dim(dev.annot)
dev.annot <- merge(dev.annot, three.agg, by="genes",all.x=T)
dim(dev.annot)
dev.annot$exon.1 <- NA
dev.annot$exon.1[grep("exon_1_", dev.annot$exons)] <- "yes"
dim(dev.annot)
head(dev.annot)

save(dev.annot, file= "~/region_bams/dev_peaks_annot/dev.peakset.annotation.Rdata")


## MAKING AN ANNOTATION FOR CG CONTENT OF REGIONS ##
dev.chrs<- (unique(seqs.dev[,1]))
dev.chrs<- gsub("KL","kl",dev.chrs)
dev.chrs<- gsub("chrUn","chrun",dev.chrs)
dev.chrs<- gsub("chrX","chrx",dev.chrs)
dev.chrs<- gsub("chrY","chry",dev.chrs)
dev.chrs<- gsub("AABR","aabr",dev.chrs)

rn6 <- list.files("/home/alussier/Rnor6/fasta", full=T)
rn6.dev.chr<- lapply(1:length(dev.chrs), function(x) rn6[grep(paste(dev.chrs[x], ".fasta",sep=""), rn6)])

rn6.fasta <- lapply(rn6.dev.chr, read.fasta)
str(rn6.fasta[[3]])

seqs.dev.cor <- gsub("KL","kl",seqs.dev)
seqs.dev.cor<- gsub("chrUn","chrun",seqs.dev.cor)
seqs.dev.cor<- gsub("chrX","chrx",seqs.dev.cor)
seqs.dev.cor<- gsub("chrY","chry",seqs.dev.cor)
seqs.dev.cor<- gsub("AABR","aabr",seqs.dev.cor)

chrs.cor <- as.factor(seqs.dev.cor[,1])
chrs.cor

sequence <- sapply(1:nrow(seqs.dev), function(x) rn6.fasta[[chrs.cor[x]]][[1]][seqs.dev[x,2]:seqs.dev[x,3]]) # pull reference sequence for feature
CpG <- sapply(1:nrow(seqs.dev), function(x) str_count(paste(sequence[[x]],collapse=''), "cg")) #total number of CpG in reference
CpG.per.bp <- sapply(1:nrow(seqs.dev), function(x) CpG[x]/length(sequence[[x]])) #CpG / base pair in feature
GC.content <- sapply(1:nrow(seqs.dev), function(x) sum(str_count(paste(sequence[[x]],collapse=''), c("c","g")))) # total number of cytosine or guanine nucleotides
GC.feature <- sapply(1:nrow(seqs.dev), function(x)GC.content[x]/length(sequence[[x]])) # GC content of feature (per base pair)

dev.CG.annot <- data.frame(seqs.dev) 
dev.CG.annot$length <- sapply(1:nrow(seqs.dev), function(x) length(sequence[[x]]))
dev.CG.annot$CpG.number <- CpG
dev.CG.annot$CpG.content <- CpG.per.bp
dev.CG.annot$GC.number <- GC.content
dev.CG.annot$GC.content <- GC.feature

dev.CG.annot <- dev.CG.annot[,c(4,1,2,3,9,5,6,7,8)]
head(dev.CG.annot)
save(dev.CG.annot, file= "~/region_bams/dev_peaks_annot/dev.peakset.CG.annotation.Rdata")

## MAKING A RAW SEQUENCES ANNOTATION ##
dev.sequences <- data.frame(seqs.dev)
dev.sequences <- dev.sequences[,c(4,1,2,3)]
dev.sequences$seq <- sapply(1:nrow(seqs.dev), function(x) paste(sequence[[x]],collapse=''))
head(dev.sequences)

save(dev.sequences, file= "~/region_bams/dev_peaks_annot/dev.peakset.sequences.Rdata")


