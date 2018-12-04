
## ANALYSIS OF STEADY-STATE GENE EXPRESSION ##
# Both days were modeled separately as follows#

#1. Surrogate Variable Analysis (SVA)
library(sva)
saline
    #S4 Vector type ExpressionSetIllumina, filtered and quantile normalized
    #contains only saline samples from 1 day (16 or 39) of one tissue (hippocampus or prefrontal cortex)
mod <- model.matrix(~prenataltreat, pData(saline))
    #primary factors - protected from sva
mod0 <- cbind(rep(1,ncol(saline)))
    #factors to correct - none in this case
svaobj.saline <- sva(saline, mod, mod0, method="irw")
    #sva will automatically calculate the number of surrogate variables

#2. Linear modeling of saline-injected animals (steady-state)
library(limma)
design <- model.matrix(~0+prenataltreat)
design <- cbind(design, svaobj.saline$sv)
    #model matrix which includes prenatal treatment and surrogate variables as cofactors
fit.saline <- lmFit(saline, design)
    #fitting data to linear model
contrast.matrix <- makeContrasts(PAE-C, PF-C, PAE-PF, levels=design)
    #setting contrasts of interest - prenatal treatments in this case
fit.saline <- contrasts.fit(fit.saline, contrast.matrix)
    #fitting model to contrasts
fit.saline <- eBayes(fit.saline)
    #compute statistics for model

#3. Identification of differentially expressed genes
top16S.PAEvC <- topTable(fit.saline, coef=1, number="Inf", p.value=1)
top16S.PFvC <- topTable(fit.saline, coef=2, number="Inf", p.value=1)
top16S.PAEvPF <- topTable(fit.saline, coef=3, number="Inf", p.value=1)
    #extract complete table of genes along with statistics corresponding to the contrast of interest
    #these tables were also used in downstream analyses - e.g. GO analysis, ErmineJ, Ingenuity Pathway Analysis
nf.saline.t25 <- topTable(fit.saline, sort.by="F", number="Inf", p.value=.25)
    #genes that are significant at F-test fdr 25% with all treatments using nestedF test at FDR <0.25%
m1 <- merge(top16S.PAEvC, top16S.PFvC, by="row.names")
m2 <- merge(m1, top16S.PAEvPF, by.x="Row.names", by.y="row.names")
m3 <- merge(nf16.t25, m2, by.x="row.names", by.y="Row.names")
pae.specific <- m3[m3$pval.EvC <0.05 & m3$pval.EvPF <0.05 & m3$pval.PFvC >0.05,]
    #genes that are differentially expressed in PAE, but not in PF vs C animals



## ANALYSIS OF GENE EXPRESSION FOLLOWING INFLAMMATORY CHALLENGE ## 
# Both days were modeled separately as follows#

#1. Surrogate Variable Analysis (SVA)
library(sva)
salineVSadjuvant 
    #S4 Vector type ExpressionSetIllumina, filtered and quantile normalized
    #contains only samples from 1 day (16 or 39) of one tissue (hippocampus or prefrontal cortex)
mod <- model.matrix(~prenataltreat*adjtreat, pData(salineVSadjuvant))
    #primary factors - protected from sva
mod0 <- cbind(rep(1,ncol(salineVSadjuvant)))
    #factors to correct - none in this case
svaobj.salineVSadjuvant <- sva(salineVSadjuvant, mod, mod0, method="irw")
    #sva will automatically calculate the number of surrogate variables


#2. Linear modeling of gene expression following saline or adjuvant injection
library(limma)
treats <- paste(salineVSadjuvant$prenataltreat, salineVSadjuvant$adjtreat, sep=".")
treats <- factor(treats, levels=c("C.S","PF.S","PAE.S","C.A","PF.A","PAE.A"))
design <- model.matrix(~0+treats)
colnames(design) <- levels(treats)
design <- cbind(design, svaobj.salineVSadjuvant$sv)
    #model matrix which includes prenatal treatment per adjuvant exposure group and surrogate variables as cofactors
 
fit.salineVSadjuvant <- lmFit(salineVSadjuvant, design)
    #fitting data to linear model
contrast.matrix <- makeContrasts(C.A-C.S, PF.A-PF.S, PAE.A-PAE.S,
                                 PAEvsC=(PAE.A-PAE.S)-(C.A-C.S),
                                 PAEvsPF=(PAE.A-PAE.S)-(PF.A-PF.S),
                                 PFvsC=(PF.A-PF.S)-(C.A-C.S),
                                 levels=design)
    #setting contrasts of interest - prenatal treatments per adjuvant exposure group in this case
fit.salineVSadjuvant <- contrasts.fit(fit.salineVSadjuvant, contrast.matrix)
    #fitting model to contrasts
fit.salineVSadjuvant <- eBayes(fit.salineVSadjuvant)
    #compute statistics for model

#3. Identification of differentially expressed genes
topC <- topTable(fit.salineVSadjuvant, coef=1, number="Inf", p.value=1)
topPF <- topTable(fit.salineVSadjuvant, coef=2, number="Inf", p.value=1)
topPAE <- topTable(fit.salineVSadjuvant, coef=3, number="Inf", p.value=1)
topPAEvC <- topTable(fit.salineVSadjuvant, coef=4, number="Inf", p.value=1)
topPAEvPF <- topTable(fit.salineVSadjuvant, coef=5, number="Inf", p.value=1)
topPFvC <- topTable(fit.salineVSadjuvant, coef=6, number="Inf", p.value=1)
    #extract complete table of genes along with statistics corresponding to the contrast of interest
    #These tables were alse used in downstream analyses - e.g. GO analysis, ErmineJ, Ingenuity Pathway Analysis
nf.salineVSadjuvant.t25 <- topTable(fit.salineVSadjuvant, sort.by="F", number="Inf", p.value=.25)
    #genes that are significant at F-test fdr 25% with all treatments using nestedF test at FDR <0.25%
m1 <- merge(top16S.PAEvC, top16S.PFvC, by="row.names")
m2 <- merge(m1, top16S.PAEvPF, by.x="Row.names", by.y="row.names")
m3 <- merge(nf.salineVSadjuvant.t25, m2, by.x="row.names", by.y="Row.names")
pae.specific <- m3[m3$pval.EvC <0.05 & m3$pval.EvPF <0.05 & m3$pval.PFvC >0.05,]
    #genes that are differentially expressed in response to adjuvant injection in PAE, but not in PF and C animals 





