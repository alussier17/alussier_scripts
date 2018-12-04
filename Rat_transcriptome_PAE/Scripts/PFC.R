

## CODE FOR DATA IN BRIEF PAPER ##

#PREFRONTAL CORTEX
library(limma)
library(beadarray)
library(gplots)
library(RColorBrewer)


setwd("/Users/Alex/Desktop/KS analysis/")
load("BSDataBsh2.RData")
PFCdesign <- read.delim("PFC_exp_parameters.txt", header=T, row.names=1)
PFC.no.out <- read.delim("PFC_exp_parameters_samplesonly_noOutliers.txt", header=T,row.names=1)
PFCdesign<- PFCdesign[order(rownames(PFCdesign)),]
head(PFCdesign)
dim(PFCdesign)

identical(colnames(exprs(BSDataBsh2)), rownames(PFCdesign))
BSData.quant <- normaliseIllumina(BSDataBsh2, method = "quantile", transform ="none")

PFC.out <- PFCdesign[!rownames(PFCdesign) %in% rownames(PFC.no.out),]

col = colorRampPalette(brewer.pal(10,"RdYlBu")[1:10])(30)

# correlation of hyb reps 
hybrep <- PFCdesign[((PFCdesign$breeding=="C41")&(PFCdesign$treatments =="C_A_16")),]
data.hybrep <- BSData.quant[,rownames(hybrep)]

phenoData(data.hybrep) <- new("AnnotatedDataFrame", data=hybrep)
identical(colnames(exprs(data.hybrep)), rownames(hybrep))

heatmap.2(cor(exprs(data.hybrep), method="pearson"), trace="none", margins=c(10,10), col=col, dendrogram="col",
          )
cor.hyb <- cor(exprs(data.hybrep), method="pearson")
quant.cor.hyb <- quantile(cor.hyb)


#correlation of amp reps
amprep <- PFCdesign[((PFCdesign$treatments=="C_S_39")&(PFCdesign$breeding=="C44")),]
data.amprep <- BSData.quant[,rownames(amprep)]

phenoData(data.amprep) <- new("AnnotatedDataFrame", data=amprep)
heatmap.2(cor(exprs(data.amprep), method="pearson"), trace="none", margins=c(15,15), cexRow = 2, cexCol = 2)
cor.amp <- cor(exprs(data.amprep), method="pearson")
quant.cor.amp <- quantile(cor.amp)


#correlation of amp and hyb reps
amp.hyb <- rbind(hybrep, amprep)
data.amp.hyb <- BSData.quant[,rownames(amp.hyb)]
phenoData(data.amp.hyb) <- new("AnnotatedDataFrame", data=amp.hyb)
amp.hyb.cor <- cor(exprs(data.amp.hyb), method="pearson")
quant.amp.hyb.cor <- quantile(amp.hyb.cor)


breaks=seq(0.90,1, by=0.0025)
col = colorRampPalette(brewer.pal(10,"RdYlBu")[1:10])(30)

amp.hyb.side <- as.matrix(c("grey", rep("black", nrow(amp.hyb)-3), "gray35","brown"))
colnames(amp.hyb.side) <- "Sample Type"

heatmap.3(cor(exprs(data.amp.hyb), method="pearson"), trace="none", margins=c(2,10), 
          cexRow = 1.5, cexCol = 2, Colv=T, Rowv=T, dendrogram="col",
          ColSideColors=amp.hyb.side,  ColSideColorsSize = 4,
          density.info="none", col=col, labCol=NA, keysize=2, KeyValueName="")
plot.new()
legend("left", legend=c("Sample hyb_rep", "Hyb_rep","Sample amp_rep", "Amp_rep"), 
       fill=c("grey","black","gray35", "brown"), border=FALSE, bty="n", y.intersp = 1, cex=1.3, x.intersp=0.35)

# correlation of all samples
phenoData(BSData.quant) <- new("AnnotatedDataFrame", data=PFCdesign)
heatmap.2(cor(exprs(BSData.quant), method="pearson"), trace="none",margins=c(7,7), dendrogram ="col", col=col)
cor.all <- (cor(exprs(BSData.quant), method="pearson"))
quant.cor.all <- quantile(cor.all)

sidecols <- as.matrix(read.delim("/Users/Alex/Desktop/DiB Paper/Treat_colours.txt", header=T))
colnames(sidecols) <- c("Sample Type","Treatment","Adjuvant","Day")
colours <-  as.matrix(read.delim("/Users/Alex/Desktop/DiB Paper/colours.txt", header=T))

source("/Users/Alex/Desktop/heatmap.3.R")
heatmap.3(cor(exprs(BSData.quant), method="pearson"), trace="none",margins=c(2,2), 
        dendrogram="col",col=col, ColSideColors=sidecols, ColSideColorsSize = 4, Rowv=T, Colv=T,
        labCol=NA, keysize=2,KeyValueName="",labRow=NA)
heatmap.3(cor(exprs(BSData.quant), method="pearson"), trace="none",margins=c(7,7), 
          col=col, ColSideColors=sidecols[], ColSideColorsSize = 4, Rowv=F, Colv=F)

plot.new()
legend("left", legend=c("Amp_Rep","Hyp_Rep","Sample","", "Con","PF","PAE","","Saline","Adjuvant","","D16","D39"), 
       fill=c("brown","black","grey", "white", "blue","green","red","white", "darkgreen","purple","white","orange","gold"), 
       border=FALSE, bty="n", y.intersp = 0.8, cex=1.3, x.intersp=0.35)



## REMOVING OUTLIERS AND REDOING NORM + CORRELATIONS

BSDataBsh <- BSDataBsh2[,rownames(PFC.no.out)]
identical(colnames(exprs(BSDataBsh)), rownames(PFC.no.out))

phenoData(BSDataBsh) <- new("AnnotatedDataFrame", data=PFC.no.out)
a <- BSDataBsh
det = calculateDetection((a), status=fData(a)$Status, negativeLabel="negative")
Detection(a) = det
expressed = apply(Detection(a) < 0.05, 1, any)
BSData.expressed = a[expressed, ]
BSData.expgenes = BSData.expressed[which(fData(BSData.expressed)$Status =="Gene"), ]

dim(exprs(BSDataBsh2)) # 23350 probes
dim(exprs(BSData.expressed)) # 20842 probes
dim(exprs(BSData.expgenes)) # 20215 probes

BSData.filt <- BSData.expgenes
BSData.filt.quant <- normaliseIllumina(BSData.filt, method = "quantile", transform ="none")

filt.cor.all <- (cor(exprs(BSData.filt.quant), method="pearson"))
quant.filt.cor.all <- quantile(filt.cor.all)

rownames(sidecols)<- rownames(PFCdesign)
sidecols.no.out <- sidecols[rownames(PFC.no.out),]

rownames(colours) <-  rownames(PFCdesign)
colours.no.out <- colours[rownames(PFC.no.out),]

heatmap.3(cor(exprs(BSData.filt.quant), method="pearson"), trace="none",margins=c(2,2), 
          col=col, Rowv=T, Colv=T, ColSideColors=sidecols.no.out[,-1], ColSideColorsSize = 4, dendrogram="col",
          labCol=NA, keysize=2,KeyValueName="",labRow=NA )
legend("left", legend=c("Con","PF","PAE","","Saline","Adjuvant","","D16","D39"), 
       fill=c("blue","green","red","white", "darkgreen","purple","white","orange","gold"), 
       border=FALSE, bty="n", y.intersp = 0.8, cex=1.3, x.intersp=0.35)
#for use with "colours no out"
legend("left", legend=c("Con","PF","PAE","","Saline","Adjuvant","","D16","D39","",
                        "Dissection 1", "Dissection 2", "Dissection 3","","Extraction 1","Extraction 2", "",
                        "Amplification 1","Amplification 2","Amplification 3","Amplification 4","","Chip 5398636011",
                        "Chip 5398636033", "Chip 5398636034","Chip 5398636035","Chip 5398636036","Chip 5398636038",
                        "Chip 5398636039","Chip 5398636044"),
       fill=c("blue","green","red","white", "darkgreen","purple","white","orange","gold","white",
              "lightblue","lightslateblue","navy","white","olivedrab1","olivedrab","white",
              "palevioletred1","palevioletred4","rosybrown1","rosybrown4","white",
              "snow1","snow3","snow4","wheat","wheat3","wheat4","yellow4","seashell4"),
       border=FALSE, bty="n", y.intersp = 1, cex=0.7, x.intersp=0.35)



# BOXPLOT OF QUANT VS NON QUANT NORM
phenoData(BSData.filt) <- new("AnnotatedDataFrame",data=PFC.no.out)
par(mfrow=c(1,2))
my.data <- BSData.filt[,order(pData(BSData.filt)$treatments)]
my.main <- "Expression for filtered PFC data"
my.names <- pData(my.data)$treatments
boxplot(exprs(my.data), main=my.main, outline=F, names=my.names, cex.axis=.6)

my.data <- BSData.filt.quant[,order(pData(BSData.filt.quant)$treatments)]
my.main <- "Expression for filtered and\nquantile normalized PFC data"
my.names <- pData(my.data)$treatments
boxplot(exprs(my.data), main=my.main, outline=F, names=my.names, cex.axis=.6)



## PCA ##
metadata <- read.delim("/Users/Alex/Desktop/KS analysis/PFC_parameters_filt2.txt", header=T, row.names=1)
meta_continuous <- metadata[,c(20,26)]
meta_continuous$RIN <- as.integer(meta_continuous$RIN)
colnames(meta_continuous) <- c("Clinical Score","RIN")
str(meta_continuous)

meta_categorical <- metadata[,c(3,5,6,7,9:12,15)]
meta_categorical$day <- as.factor(meta_categorical$day)
meta_categorical$set <- as.factor(meta_categorical$set)
meta_categorical$dissection <- as.factor(meta_categorical$dissection)
meta_categorical$extraction <- as.factor(meta_categorical$extraction)
meta_categorical$amplification <- as.factor(meta_categorical$amplification)
meta_categorical$chipbarcode <- as.factor(meta_categorical$chipbarcode)
str(meta_categorical)
colnames(meta_categorical) <- c("Prenatal treatment","Adjuvant","Day","Set","Dissection","Extraction",
                                "Amplification","Chip","Developed AA")
dim(meta_categorical)

PCA_full<- princomp(exprs(BSData.filt.quant), cor=F)
Loadings <- as.data.frame(unclass(PCA_full$loadings))
vars <- PCA_full$sdev^2
Importance<-vars/sum(vars)
heat_scree_plot(Loadings, Importance)


## DATA RESTRUCTURING ##
design2 <- read.delim("PFC_parameters_filt2.txt", header=T, row.names=1)

BSData <- new("ExpressionSet", phenoData = as(design2, "AnnotatedDataFrame"), exprs = as.matrix(exprs(BSData.filt.quant)))

dim(BSData)
# 85 samples, 20215 probes
str(pData(BSData))

BSData$DevelopedAA <- factor(BSData$DevelopedAA, c("Saline","No","Yes"))
BSData$AAcategory <- factor(BSData$AAcategory, c("Saline","None","Mild","Moderate"))
BSData$adjtreat <- relevel(BSData$adjtreat, "S")
BSData$prenataltreat <- factor(BSData$prenataltreat, c("C","PF","E"))
BSData$amplification <- factor(BSData$amplification)
BSData$day <- factor(BSData$day, c("16","39"))
BSData$treatDev1 <- factor(BSData$treatDev1, c("Saline_C","Saline_PF","Saline_E","No_C","No_PF","No_E","Yes_C","Yes_PF","Yes_E"))
BSData$treatDev2 <- factor(BSData$treatDev2, c("C_Saline","C_No","C_Yes","PF_Saline","PF_No","PF_Yes","E_Saline","E_No","E_Yes"))
BSData$dissection <- factor(BSData$treatDev2, c("1","2","3"))

## SVA ## 
library(beadarray)
library(sva)
library(limma)

dat16 <- BSData[,BSData$day=="16"]
dat39 <- BSData[,BSData$day=="39"]

datS <- BSData[,BSData$adjtreat=="S"]
datS.16 <- datS[,datS$day=="16"]
datS.39 <- datS[,datS$day=="39"]

#steady-state day 16
svadat.S16 <- exprs(datS.16)
mod.S16 <- model.matrix(~prenataltreat, pData(datS.16))
dim(mod.S16)
mod0.S16 <- cbind(rep(1,ncol(exprs(datS.16))))
svaobj.S16 <- sva(svadat.S16, mod.S16, mod0.S16, method="irw")
dim(svaobj.S16$sv) #6 sv
#steady-state day 39
svadat.S39 <- exprs(datS.39)
mod.S39 <- model.matrix(~prenataltreat, pData(datS.39))
dim(mod.S39)
mod0.S39 <- cbind(rep(1,ncol(exprs(datS.39))))
svaobj.S39 <- sva(svadat.S39, mod.S39, mod0.S39, method="irw")
dim(svaobj.S39$sv) #4 sv
#steady state
svadat.S <- exprs(datS)
mod.S <- model.matrix(~prenataltreat, pData(datS))
dim(mod.S)
mod0.S <- cbind(rep(1,ncol(exprs(datS)))) # make empyt matrix with same number of rows as mod
svaobj.S <- sva(svadat.S, mod.S, mod0.S, method="irw")
dim(svaobj.S$sv) #10 sv

#Saline VS Adjuvant - Day 16
svadat.SvA.16 <- exprs(dat16)
mod.SvA.16 <- model.matrix(~prenataltreat*adjtreat, pData(dat16))
dim(mod.SvA.16)
mod0.SvA.16 <- cbind(rep(1,ncol(exprs(dat16))))
svaobj.SvA.16 <- sva(svadat.SvA.16, mod.SvA.16, mod0.SvA.16, method="irw")
dim(svaobj.SvA.16$sv) #15 sv
#Saline VS Adjuvant - Day 39
svadat.SvA.39 <- exprs(dat39)
mod.SvA.39 <- model.matrix(~prenataltreat*adjtreat, pData(dat39))
dim(mod.SvA.39)
mod0.SvA.39 <- cbind(rep(1,ncol(exprs(dat39))))
svaobj.SvA.39 <- sva(svadat.SvA.39, mod.SvA.39, mod0.SvA.39, method="irw")
dim(svaobj.SvA.39$sv) #12 sv
# entire PFC dataset
mod.all <- model.matrix(~prenataltreat, pData(BSData))
dim(mod.all)
mod0.all <- cbind(rep(1,ncol(exprs(BSData))))
sva.all <- sva(exprs(BSData), mod.all, mod0.all, method="irw")
t <- sva.all$sv


# SALINE HITS # 
eset <- dat16[,dat16$adjtreat=="S"]
treats <- paste(eset$prenataltreat)
treats <- factor(treats, levels=c("C","PF","E"))
design <- model.matrix(~0+treats)
colnames(design) <- levels(treats)
sv <- svaobj.S16$sv
dim(sv)
colnames(sv) <- c("sv1","sv2","sv3","sv4","sv5","sv6")
design <- cbind(design, sv)

fit <- lmFit(eset, design)
contrast.matrix <- makeContrasts(E-C, PF-C, E-PF, levels=design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
fit16.pntS <- fit2

results16.pntS <- decideTests(fit16.pntS, adjust.method="BH", p.value=.05)
summary(results16.pntS)
results16.pntS.20 <- decideTests(fit16.pntS, adjust.method="BH", p.value=.2)
summary(results16.pntS.20)
results16.pntS.25 <- decideTests(fit16.pntS, adjust.method="BH", p.value=.25)
summary(results16.pntS.25)
results16.pntS.50 <- decideTests(fit16.pntS, adjust.method="BH", p.value=.5)
summary(results16.pntS.50)
top16S.EvC <- topTable(fit16.pntS, coef=1, number="Inf", p.value=1)
top16S.PFvC <- topTable(fit16.pntS, coef=2, number="Inf", p.value=1)
top16S.EvPF <- topTable(fit16.pntS, coef=3, number="Inf", p.value=1)

par(mfrow=c(1,1))
par(lwd=3)
plot(0, xlab="P-value", ylab="Average density", type="n", xlim=c(0,1), ylim=c(0,2), 
     main="P-value distributions for D16 PFC\nPrenatal diet comparisons (Saline SVA)", xaxs="i", yaxs="i")
lines(density(top16S.EvC$P.Value), col=cols[1], lty=1)
lines(density(top16S.EvPF$P.Value), col=cols[2], lty=2)
lines(density(top16S.PFvC$P.Value), col=cols[3], lty=6)
legend("topright",c("PAE-C","PAE-PF","PF-C"), col=c("red","orange","green"), lty=c(1,2,6))
venn <- results16.pntS.20
colnames(venn) <- c("E-C","PF-C","E-PF")
vennDiagram(venn)
mtext("PFC D16: Genes DE among prenatal diet groups,\nat FDR of 20%", line=1, cex=1.3)



## NOT WORKING AT ALL ##
library(Rarity)
library(lattice)
library(car)
test <-cbind(sva.all$sv, meta_categorical$dissection, meta_categorical$extraction,meta_categorical$amplification,
             meta_categorical$chipbarcode)
pairs(test)

Loadings <- as.data.frame(unclass(sva.all$sv))
vars <- (sapply(1:ncol(sva.all$sv), function(x) sd(sva.all$sv[,x])))^2
Importance<-vars/sum(vars)
heat_scree_plot(Loadings, Importance)

cor(~sva.all$sv + meta_categorical)


sapply(1:ncol(sva.all$sv), function(x) sd(sva.all$sv[,x])) # sd(sva.all$sv[,x]))
?apply

sva.all$sv[,1]

sd(sva.all$sv[,15])
dim(sva.all$sv)



