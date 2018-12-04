
## CODE FOR DATA IN BRIEF PAPER ##

#HIPPOCAMPUS
library(limma)
library(beadarray)
library(gplots)
library(RColorBrewer)

setwd("/Users/Alex/Desktop/KS analysis/")
load("BSDataBsh.hipp.RData")
dim(exprs(BSDataBsh.hipp))
design <- read.delim("Hipp_parameters_all.txt", header=T, row.names=1)
dim(design)

identical(colnames(exprs(BSDataBsh.hipp)), rownames(design))
BSData.hipp.quant <- normaliseIllumina(BSDataBsh.hipp, method = "quantile", transform ="none")

col = colorRampPalette(brewer.pal(10,"RdYlBu")[1:10])(30)
colours <- as.matrix(read.delim("/Users/Alex/Desktop/DiB Paper/colours_hipp.txt"))
source("/Users/Alex/Desktop/heatmap.3.R")

# correlation of hyb reps 
hybrep <- design[(design$sampletype=="hybrep"),]
hybrep <- rbind(hybrep, design[(design$breeding=="PF30")&(design$sampletype=="sample")&(design$treatments=="PF_A_39"),],
                design[(design$breeding=="PF28")&(design$sampletype=="sample")&(design$treatments=="PF_S_16"),],
                design[(design$breeding=="E11")&(design$sampletype=="sample")&(design$treatments=="E_A_39"),],
                design[(design$breeding=="C46")&(design$sampletype=="sample")&(design$treatments=="C_A_39"),])

data.hybrep <- BSData.hipp.quant[,rownames(hybrep)]
phenoData(data.hybrep) <- new("AnnotatedDataFrame", data=hybrep)
identical(colnames(exprs(data.hybrep)), rownames(hybrep))

test <- cbind(rownames(design), colours)
test.2 <- merge(test, hybrep, by.x="V1", by.y="row.names")
cols.hybrep <- as.matrix(cbind(test.2[,2:5], 
                c("darkgreen","darkgreen","olivedrab","olivedrab","darkred","darkblue","darkgreen","darkblue","darkred",
                "olivedrab","darkred","darkblue","darkgreen","olivedrab","darkred","darkblue")))
colnames(cols.hybrep) <- c(colnames(cols.hybrep[,1:4]), "Replicate Group")


heatmap.3(cor(exprs(data.hybrep), method="pearson"), trace="none", margins=c(10,10), col=col, dendrogram="col",
          ColSideColors=cols.hybrep[,c(1,5)],  ColSideColorsSize = 3, Rowv=T, Colv=T)
legend("left", legend=c("Hyp_Rep","Sample","","Group 1",
                        "Group 2", "Group 3", "Group 4"), 
              fill=c("black","grey", "white", "darkgreen","olivedrab","darkred","darkblue"), 
              border=FALSE, bty="n", y.intersp = 1, cex=0.7, x.intersp=0.35)

cor.hyb <- cor(exprs(data.hybrep), method="pearson")
quant.cor.hyb <- quantile(cor.hyb)
quant.cor.hyb
quant.hyb.1 <- quantile(cor.hyb[c(1,2,7,13),c(1,2,7,13)])
quant.hyb.2 <- quantile(cor.hyb[c(3,4,10,14),c(3,4,10,14)])
quant.hyb.3 <- quantile(cor.hyb[c(5,9,11,15),c(5,9,11,15)])
quant.hyb.4 <- quantile(cor.hyb[c(6,8,12,16),c(6,8,12,16)])

hybrep.ordered <- rbind(hybrep[c(1,2,7,13),], hybrep[c(3,4,10,14),], hybrep[c(5,9,11,15),], hybrep[c(6,8,12,16),])
hybrep.cols.ordered <- rbind(cols.hybrep[c(1,2,7,13),], cols.hybrep[c(3,4,10,14),], 
                             cols.hybrep[c(5,9,11,15),], cols.hybrep[c(6,8,12,16),])
colnames(hybrep.cols.ordered) <- c("Sample Type", "Treatment","Adj","Day","Replicate Group")


hybrep.data.ordered <- BSData.hipp.quant[,rownames(hybrep.ordered)]
phenoData(hybrep.data.ordered ) <- new("AnnotatedDataFrame", data=hybrep.ordered)
identical(colnames(exprs(hybrep.data.ordered)), rownames(hybrep.ordered))

par(cex.axis = 2)
heatmap.3(cor(exprs(hybrep.data.ordered ), method="pearson"), trace="none", margins=c(2,10), 
          col=col, dendrogram="col", ColSideColors=hybrep.cols.ordered[,c(1,5)],  
          ColSideColorsSize = 4, Rowv=T, Colv=T, keysize=2, labCol=NA, KeyValueName="", cexRow=1.5)
plot.new()
legend("bottomleft", legend=c("Hyb_rep","Sample","","Group 1",
                        "Group 2", "Group 3", "Group 4"), 
       fill=c("black","grey", "white","darkgreen","olivedrab","darkred","darkblue"), 
       border=FALSE, bty="n", y.intersp = 0.8, cex=1.3, x.intersp=0.35)



# correlation of all samples
phenoData(BSData.hipp.quant) <- new("AnnotatedDataFrame", data=design)
heatmap.2(cor(exprs(BSData.hipp.quant), method="pearson"), trace="none",margins=c(7,7), dendrogram ="col", col=col)
cor.all <- (cor(exprs(BSData.hipp.quant), method="pearson"))
quant.cor.all <- quantile(cor.all)
quant.cor.all
par(cex.axis = 2)
colnames(colours) <- c("Sample Type","Treatment","Adjuvant","Day")

par(cex.axis = 2)
heatmap.3(cor(exprs(BSData.hipp.quant), method="pearson"), trace="none",margins=c(2,2), 
          dendrogram="col",col=col, ColSideColors=colours, ColSideColorsSize = 4, 
          Rowv=T, Colv=T, labCol=NA, keysize=2,KeyValueName="",labRow=NA)
plot.new()
legend("left", legend=c("Hyp_rep","Sample","", "Con","PF","PAE",
                        "","Saline","Adjuvant","","D16","D39"), 
       fill=c("black","grey", "white", "blue","green","red","white", 
              "darkgreen","purple","white","orange","gold"), 
       border=FALSE, bty="n", y.intersp = 0.8, cex=1.3, x.intersp=0.35)


## REMOVING HYBRIDIZATION REPLICATES ## 
design.samples <- design[design$sampletype=="sample",]
dim(design.samples)
BSDataBsh <- BSDataBsh.hipp[,rownames(design.samples)]
dim(exprs(BSDataBsh)) ## now 84 columns, due to filtering
identical(colnames(exprs(BSDataBsh)), rownames(design.samples))
phenoData(BSDataBsh) <- new("AnnotatedDataFrame", data=design.samples)

a <- BSDataBsh
det = calculateDetection((a), status=fData(a)$Status, negativeLabel="negative")
Detection(a) = det
expressed = apply(Detection(a) < 0.05, 1, any)
BSData.expressed = a[expressed, ]
BSData.expgenes = BSData.expressed[which(fData(BSData.expressed)$Status =="Gene"), ]

dim(exprs(BSDataBsh)) # 23350 probes
dim(exprs(BSData.expressed)) # 20673 probes
dim(exprs(BSData.expgenes)) # 20069 probes

BSData.filt <- BSData.expgenes
BSData <- normaliseIllumina(BSData.filt, method = "quantile", transform ="none")

filt.cor.all <- (cor(exprs(BSData), method="pearson"))
quant.filt.cor.all <- quantile(filt.cor.all)

rownames(colours) <- rownames(design)
cols.no.out <- colours[rownames(design.samples),]

heatmap.3(cor(exprs(BSData), method="pearson"), trace="none",margins=c(2,2), 
          col=col, Rowv=T, Colv=T, ColSideColors=cols.no.out[,-1], ColSideColorsSize = 4, dendrogram="col",
          labCol=NA, keysize=2,KeyValueName="",labRow=NA)
heatmap.3(cor(exprs(BSData), method="pearson"), trace="none",margins=c(7,7), 
          col=col, Rowv=F, Colv=F, ColSideColors=cols.no.out, ColSideColorsSize = 3, dendrogram="col")
plot.new()
legend("left", legend=c("Con","PF","PAE","","Saline","Adjuvant","","D16","D39"), 
       fill=c("blue","green","red","white", "darkgreen","purple","white","orange","gold"), 
       border=FALSE, bty="n", y.intersp = 0.8, cex=1.3, x.intersp=0.35)

detail_col <- as.matrix(read.delim("/Users/Alex/Desktop/DiB Paper/hipp_detail_colours.txt", header=T, row.names=1))
detail_col.samples <- detail_col[rownames(design.samples),]

heatmap.3(cor(exprs(BSData), method="pearson"), trace="none",margins=c(7,7), 
          col=col, Rowv=T, Colv=T, ColSideColors=detail_col.samples, ColSideColorsSize = 4, dendrogram="col")
heatmap.3(cor(exprs(BSData), method="pearson"), trace="none",margins=c(7,7), 
          col=col, Rowv=F, Colv=F, ColSideColors=detail_col.samples, ColSideColorsSize = 4, dendrogram="col")
plot.new()
par(mar=c(1,1,1,1))
legend("left", legend=c("Con","PF","PAE","","Saline","Adjuvant","","D16","D39","",
                        "Dissection 1", "Dissection 2", "Dissection 3","","Extraction 1","Extraction 2", "",
                        "Amplification 1","Amplification 2","Amplification 3","Amplification 4","","Chip 5398636045",
                        "Chip 5398636049", "Chip 5398636050","Chip 5398636051","Chip 5398636052","Chip 5398636055",
                        "Chip 5430261013","Chip 5430261014"),
       fill=c("blue","green","red","white", "darkgreen","purple","white","orange","gold","white",
              "lightblue","lightslateblue","navy","white","olivedrab1","olivedrab","white",
              "palevioletred1","palevioletred4","rosybrown1","rosybrown4","white",
              "snow1","snow3","snow4","wheat","wheat3","wheat4","yellow4","seashell4"),
       border=FALSE, bty="n", y.intersp = 1, cex=0.7, x.intersp=0.35)




# BOXPLOT OF QUANT VS NON QUANT NORM
par(mfrow=c(1,2))
my.data <- BSData.filt[,order(pData(BSData.filt)$treatments)]
my.main <- "Expression for filtered HPC data"
my.names <- pData(my.data)$treatments
boxplot(exprs(my.data), main=my.main, outline=F, names=my.names, cex.axis=.6)

my.data <- BSData[,order(pData(BSData)$treatments)]
my.main <- "Expression for filtered and\nquantile normalized HPC data"
my.names <- pData(my.data)$treatments
boxplot(exprs(my.data), main=my.main, outline=F, names=my.names, cex.axis=.6)
par(mfrow=c(1,1))

## PCA ##
metadata <- read.delim("/Users/Alex/Desktop/KS analysis/RIN comparison/Hipp_allparameters_samplesonly.txt", header=T, row.names=1)
colnames(metadata)
meta_continuous <- metadata[,c(30,17)]
meta_continuous$RIN <- as.integer(meta_continuous$RIN)
colnames(meta_continuous) <- c("Clinical Score","RIN")
str(meta_continuous)

meta_categorical <- metadata[,c(8,9,10,11,13:15,2,29)]
meta_categorical$Day <- as.factor(meta_categorical$Day)
meta_categorical$Set <- as.factor(meta_categorical$Set)
meta_categorical$Dissection <- as.factor(meta_categorical$Dissection)
meta_categorical$Extraction <- as.factor(meta_categorical$Extraction)
meta_categorical$Amplification <- as.factor(meta_categorical$Amplification)
meta_categorical$Array <- as.factor(meta_categorical$Array)
colnames(meta_categorical) <- c("Prenatal treatment","Adjuvant","Day","Set","Dissection","Extraction",
                                "Amplification","Chip","Developed AA")
str(meta_categorical)
dim(meta_categorical)

PCA_full<- princomp(exprs(BSData), cor=F)
Loadings <- as.data.frame(unclass(PCA_full$loadings))
vars <- PCA_full$sdev^2
Importance<-vars/sum(vars)
heat_scree_plot(Loadings, Importance)

eset <- BSData
mvals <- is.na(exprs(eset))
sum(mvals) ## no missing values in my data

pcs <- prcomp(exprs(eset), center=T,scale=T)
x <- summary(pcs)
str(x)
screeplot(pcs)
plot(x$importance[2,], ylab="Proportion of Variance", xlab="Principal Component", main="Proportion of Variance attributable to each PC")
pData(eset) <- cbind(pData(eset), pcs$rotation[sampleNames(eset),1:20])
BSData.pcs <- eset
plot(pData(eset)[,c(5:8, 16:30)], main="Principal Components vs Prenatal treat")
colnames(pData(eset))

#Limma
eset$adjtreat <- relevel(eset$adjtreat, "S")
eset$prenataltreat <- factor(eset$prenataltreat, c("C","PF","E"))
eset$amplification <- factor(eset$amplification)
eset$chipnumber <- factor(eset$chipnumber)
eset$dissection <- factor(eset$dissection)

# SVA # 


