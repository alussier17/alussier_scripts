#creating read-in file for the summary stats format

# set working directory
setwd("YOURPATH")

#-------------------------------------------------------------------------------
# MENTAL HEALTH OUTCOMES
#-------------------------------------------------------------------------------

### ADHD -----------------------------------------------------------------------

temp <- read.table("ADHD2022_iPSYCH_deCODE_PGC.meta.gz", 
                   header=T, nrows=2)
readLines("ADHD2022_iPSYCH_deCODE_PGC.meta.gz",n=2)
temp

adhd2 <- data.frame(psych = "ADHD2022",
                    filename= "ADHD2022_iPSYCH_deCODE_PGC.meta.gz",
                    sep =" ",
                    snp_col = "SNP",
                    beta_col = "OR",
                    se_col = "SE",
                    effect_allele_col="A1",
                    other_allele_col = "A2",
                    eaf_col = "FRQ_A_38691",
                    pval_col = "P",
                    note = "OR")

### Anorexia nervosa -----------------------------------------------------------

temp <- read.table("pgcAN2.2019-07.vcf.tsv.gz", 
                   header=T, nrows=2)
readLines("pgcAN2.2019-07.vcf.tsv.gz", n=2)
temp

anorexia <- data.frame(psych = "anorexia",
                       filename= "pgcAN2.2019-07.vcf.tsv.gz",
                       sep ="\t",
                       snp_col = "ID",
                       beta_col = "BETA",
                       se_col = "SE",
                       effect_allele_col="REF",
                       other_allele_col = "ALT",
                       eaf_col = "",
                       pval_col = "PVAL",
                       note = "beta")

### Anxiety disorders ----------------------------------------------------------

temp <- read.table("anxiety.meta.full.cc.tbl.gz", 
                   header=T, nrows=2)
temp
readLines("anxiety.meta.full.cc.tbl.gz",n=1)

anxiety <- data.frame(psych = "anxiety",
                      filename= "anxiety.meta.full.cc.tbl.gz",
                      sep ="\t",
                      snp_col = "SNPID",
                      beta_col = "Effect",
                      se_col = "StdErr",
                      effect_allele_col="Allele1",
                      other_allele_col = "Allele2",
                      eaf_col = "Freq1",
                      pval_col = "P.value",
                      note = "")

### Autism spectrum disorder ---------------------------------------------------

temp <- read.table("iPSYCH-PGC_ASD_Nov2017.gz", 
                   header=T, nrows=2)
temp
readLines("iPSYCH-PGC_ASD_Nov2017.gz",n=1)

asd <- data.frame(psych = "asd",
                  filename= "iPSYCH-PGC_ASD_Nov2017.gz",
                  sep="\t",
                  snp_col = "SNP",
                  beta_col = "OR",
                  se_col = "SE",
                  effect_allele_col="A1",
                  other_allele_col = "A2",
                  eaf_col = "",
                  pval_col = "P",
                  note = "OR")

### Bipolar disorder -----------------------------------------------------------

temp <- read.table("pgc-bip2021-all.vcf.tsv.gz", 
                   header=T, nrows=2)
temp
readLines("pgc-bip2021-all.vcf.tsv.gz",n=1)

bipolar <-  data.frame(psych = "bipolar",
                       filename= "pgc-bip2021-all.vcf.tsv.gz",
                       sep="\t",
                       snp_col = "ID",
                       beta_col = "BETA",
                       se_col = "SE",
                       effect_allele_col="A1",
                       other_allele_col = "A2",
                       eaf_col = "",
                       pval_col = "PVAL",
                       note = "beta")

### Cannabis use disorder ------------------------------------------------------

temp <- data.table::fread("EUR (main)/CUD_EUR_full_public_11.14.2020.gz", 
                          header=T, nrow=2)

# temp$Neff <- 4/(1/temp$N_CAS + 1/temp$N_CON)
# temp$sigma <- 1/(1 + (temp$Z * temp$Z)/temp$Neff)
# temp$beta <- temp$Z*sqrt(temp$sigma/temp$Neff)
# temp$se <- sqrt(temp$sigma/temp$Neff)
# temp[temp$P == min(temp$P),]
# write.table(temp, file = "EUR (main)/CUD_EUR_full_public_11.14.2020_recoded.txt",
#             sep ="\t", quote=F, row.names=F, col.names=T)
# 
# 
CUD <- data.frame(pheno = "cannabis",
                  filename= "CUD_EUR_full_public_11.14.2020_recoded.txt",
                  sep ="\t",
                  snp_col = "SNP",
                  beta_col = "beta",
                  se_col = "se",
                  effect_allele_col="A1",
                  other_allele_col = "A2",
                  eaf_col = "",
                  pval_col = "P",
                  note = "recoded")

### Depression -----------------------------------------------------------------

temp <- read.table("daner_MDDwoBP_20201001_2015iR15iex_HRC_MDDwoBP_iPSYCH2015i_UKBtransformed_Wray_FinnGen_MVPaf_HRC_MAF01.gz", 
                   header=T, nrows=2)
readLines("daner_MDDwoBP_20201001_2015iR15iex_HRC_MDDwoBP_iPSYCH2015i_UKBtransformed_Wray_FinnGen_MVPaf_HRC_MAF01.gz",n=2)
temp

depression2023 <- data.frame(psych = "depression2023",
                             filename= "daner_MDDwoBP_20201001_2015iR15iex_HRC_MDDwoBP_iPSYCH2015i_UKBtransformed_Wray_FinnGen_MVPaf_HRC_MAF01.gz",
                             sep =" ",
                             snp_col = "SNP",
                             beta_col = "OR",
                             se_col = "SE",
                             effect_allele_col="A1",
                             other_allele_col = "A2",
                             eaf_col = "FRQ_A_294322",
                             pval_col = "P",
                             note = "OR")

### Obsessive-compulsive disorder ----------------------------------------------

temp <- read.table("ocd_aug2017.gz", 
                   header=T, nrows=2)
readLines("ocd_aug2017.gz",n=2)
temp

ocd <- data.frame(psych = "ocd",
           filename= "ocd_aug2017.gz",
           sep="\t",
           snp_col = "SNP",
           beta_col = "OR",
           se_col = "SE",
           effect_allele_col="A1",
           other_allele_col = "A2",
           eaf_col = "",
           pval_col = "P",
           note = "")

### Opioid dependence ----------------------------------------------------------

temp <- data.table::fread("opi.DEPvEXP_EUR.noAF.tbl",
                          header=F, skip =1)
colnames(temp) <- c("rsID","Allele1","Allele2","Weight","Zscore","P-value","HetIsq",
                    "HetChiSq","HetDf", "HetPVal","Total_N","Total_NCase","TotalNControl", "NA")
head(temp)
#weight is the effective N
temp$sigma <- 1/(1 + (temp$Zscore * temp$Zscore)/temp$Weight)
temp$beta <- temp$Zscore*sqrt(temp$sigma/temp$Weight)
temp$se <- sqrt(temp$sigma/temp$Weight)
head(temp)
write.table(temp, file = "opi.DEPvEXP_EUR.noAF_recoded.txt",
            sep ="\t", quote=F, row.names=F, col.names=T)

opioidDepVexp <- data.frame(psych = "opioidDepVexp",
                            filename= "opi.DEPvEXP_EUR.noAF_recoded.txt",
                            sep ="\t",
                            snp_col = "rsID",
                            beta_col = "beta",
                            se_col = "se",
                            effect_allele_col="Allele1",
                            other_allele_col = "Allele2",
                            eaf_col = "",
                            pval_col = "P-value",
                            note = "recoded")


load("psych_stat_summary_2023-05-02.Rdata")
psych <- psych[!psych$psych =="CAD2022",]
psych <- rbind(psych, cad2022, opioidDepVexp)
save(psych, file= "psych_stat_summary_2023-05-08.Rdata")

### PTSD -----------------------------------------------------------------------

temp <- read.table("pts_eur_freeze2_overall.results.gz", 
                   header=T, nrows=2)
readLines("pts_eur_freeze2_overall.results.gz",n=2)
temp

ptsd <-  data.frame(psych = "ptsd",
                    filename= "pts_eur_freeze2_overall.results.gz",
                    sep="\t",
                    snp_col = "SNP",
                    beta_col = "OR",
                    se_col = "SE",
                    effect_allele_col="A1",
                    other_allele_col = "A2",
                    eaf_col = "FRQ_A_23212",
                    pval_col = "P",
                    note = "OR")

### Schizophrenia --------------------------------------------------------------

temp <- read.table("PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz", 
                   header=T, nrows=2)
readLines("PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz",n=2)
temp

scz <- data.frame(psych = "scz",
                  filename= "PGC3_SCZ_wave3.european.autosome.public.v3.vcf.tsv.gz",
                  sep="\t",
                  snp_col = "ID",
                  beta_col = "BETA",
                  se_col = "SE",
                  effect_allele_col="A1",
                  other_allele_col = "A2",
                  eaf_col = "",
                  pval_col = "PVAL",
                  note = "beta")

### Suicide attempt ------------------------------------------------------------

temp <- read.table("EUR (main)/daner_model2_062620_eur.neff.qc2.80.gz", 
                   header=T, nrows=2)
readLines("EUR (main)/daner_model2_062620_eur.neff.qc2.80.gz",n=2)
temp

suicide2020 <- data.frame(psych = "suicide2020",
                          filename= "daner_model2_062620_eur.neff.qc2.80.gz",
                          sep ="\t",
                          snp_col = "SNP",
                          beta_col = "OR",
                          se_col = "SE",
                          effect_allele_col="A1.y",
                          other_allele_col = "A2.y",
                          eaf_col = "FRQ_A_26590",
                          pval_col = "P",
                          note = "OR")


psych <- rbind(psych, suicide2020)
save(psych, file = "EUR (main)/main_outcomes_2023-05-12.Rdata")

### Tourette syndrome ----------------------------------------------------------

temp <- read.table("TS_Oct2018.gz", 
                   header=T, nrows=2)
readLines("TS_Oct2018.gz",n=2)
temp

tourette <- data.frame(psych = "tourette",
                       filename= "TS_Oct2018.gz",
                       sep=" ",
                       snp_col = "SNP",
                       beta_col = "OR",
                       se_col = "SE",
                       effect_allele_col="A1",
                       other_allele_col = "A2",
                       eaf_col = "",
                       pval_col = "P",
                       note = "OR")

#-------------------------------------------------------------------------------
### PHYSICAL HEALTH OUTCOMES
#-------------------------------------------------------------------------------

### Asthma ---------------------------------------------------------------------

temp <- read.table("EUR (main)/Asthma_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz", 
                   header=F, nrows=2)
colnames(temp) <- limma::strsplit2(readLines("EUR (main)/Asthma_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz", n=1),
                                   "\t")
temp
readLines("EUR (main)/Asthma_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz",n=1)

asthma <- data.frame(pheno = "asthma",
                     filename= "Asthma_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz",
                     sep ="\t",
                     snp_col = "rsid",
                     beta_col = "inv_var_meta_beta",
                     se_col = "inv_var_meta_sebeta",
                     effect_allele_col="REF",
                     other_allele_col = "ALT",
                     eaf_col = "all_meta_AF",
                     pval_col = "inv_var_meta_p",
                     note = "beta")

### Chronic kidney disease -----------------------------------------------------

temp <- read.table("EUR (main)/CKD_overall_EA_JW_20180223_nstud23.dbgap.txt.gz", 
                   header=T, nrows=2)
temp
readLines("EUR (main)/CKD_overall_EA_JW_20180223_nstud23.dbgap.txt.gz",n=1)

CKD <- data.frame(pheno = "CKD",
                  filename= "CKD_overall_EA_JW_20180223_nstud23.dbgap.txt.gz",
                  sep =" ",
                  snp_col = "RSID",
                  beta_col = "Effect",
                  se_col = "StdErr",
                  effect_allele_col="Allele1",
                  other_allele_col = "Allele2",
                  eaf_col = "Freq1",
                  pval_col = "P.value",
                  note = "beta")

### COPD -----------------------------------------------------------------------

temp <- read.table("EUR (main)/COPD_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz", 
                   header=F, nrows=2)
colnames(temp) <- limma::strsplit2(readLines("EUR (main)/COPD_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz", n=1),
                                   "\t")
temp

readLines("EUR (main)/CUD_EUR_full_public_11.14.2020.gz",n=1)

COPD <- data.frame(pheno = "COPD",
                   filename= "COPD_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz",
                   sep ="\t",
                   snp_col = "rsid",
                   beta_col = "inv_var_meta_beta",
                   se_col = "inv_var_meta_sebeta",
                   effect_allele_col="REF",
                   other_allele_col = "ALT",
                   eaf_col = "all_meta_AF",
                   pval_col = "inv_var_meta_p",
                   note = "beta")

### Coronary artery disease ----------------------------------------------------

# head(temp)
# temp$MarkerName2 <- strsplit2(temp$MarkerName, "_")[,1]
# temp16 <- temp[temp$CHR==16,]
# temp16$MarkerName2 <- limma::strsplit2(temp16$MarkerName, "_")[,1]
# head(temp16)


#mQTLs 
load("YOURPATH/exposureSNPs_goDMC_2023-05-02A.Rdata")
dnam <- exposureSNPs
load("YOURPATH/exposureSNPs_goDMC_2023-05-02B.Rdata")
dnam <- rbind(dnam, exposureSNPs)
load("YOURPATH/exposureSNPs_goDMC_2023-05-02B.Rdata")
dnam <- rbind(dnam, exposureSNPs)
load("YOURPATH/exposureSNPs_goDMC_2023-05-10A.Rdata")
dnam <- rbind(dnam, exposureSNPs)


dnam$MarkerName <- paste(gsub("chr","", limma::strsplit2(dnam$snp, ":")[,1]),
                         limma::strsplit2(dnam$snp, ":")[,2], sep=":")
dnam$CHR <- limma::strsplit2(dnam$MarkerName, ":")[,1]
head(dnam)
unique(dnam$CHR)

for(i in unique(dnam$CHR)){
  if(i == unique(dnam$CHR)[1]){cadSubset <- data.frame()}
  print(i)
  tempChr <- temp[temp$CHR == i,]
  tempChr$MarkerName2 <- limma::strsplit2(tempChr$MarkerName, "_")[,1]
  tempChr <- tempChr[tempChr$MarkerName2 %in% dnam$MarkerName,]
  tempChr$rsid <- dnam$SNP[match(tempChr$MarkerName2, dnam$MarkerName)]
  cadSubset <- rbind(cadSubset, tempChr)
}
dim(cadSubset)
head(cadSubset)
write.table(cadSubset, file = "CAD_GWAS_primary_discovery_meta_subset.txt", sep ="\t",
            quote=F, row.names=F, col.names=T)

cad2022 <- data.frame(psych = "CAD2022",
                      filename= "CAD_GWAS_primary_discovery_meta_subset.txt",
                      sep ="\t",
                      snp_col = "rsid",
                      beta_col = "Effect",
                      se_col = "StdErr",
                      effect_allele_col="Allele1",
                      other_allele_col = "Allele2",
                      eaf_col = "Freq1",
                      pval_col = "P.value",
                      note = "recoded SNPs")


### Obesity --------------------------------------------------------------------

#Obesity - fixed firt lines
temp <- read.table("EUR (main)/Obesity.vcf.gz", 
                   header=F, nrows=2)
colnames(temp) <- limma::strsplit2(readLines("EUR (main)/Obesity.vcf.gz", n=110)[110],
                                   "\t")
temp

# obesity <- data.table::fread("EUR (main)/Obesity.vcf.gz", skip =109) 
# head(obesity[,])
# obesity2 <- cbind(obesity[,1:8],
#                   limma::strsplit2(obesity$`UKB-b-15541`,":"))
# colnames(obesity2)[9:13] <- c(limma::strsplit2(obesity$FORMAT[1],":"))
# obesity2$pval <- 10^-as.numeric(obesity2$LP)
# obesity2 <- obesity2[,-13]
# write.table(obesity2, file = "EUR (main)/Obesity_recoded.txt",sep="\t",
#             col.names=T, row.names = F, quote=F)
# rm(obesity, obesity2)

obesity <- data.frame(pheno = "obesity",
                      filename= "Obesity_recoded.txt",
                      sep ="\t",
                      snp_col = "ID",
                      beta_col = "ES",
                      se_col = "SE",
                      effect_allele_col="REF",
                      other_allele_col = "ALT",
                      eaf_col = "AF",
                      pval_col = "pval",
                      note = "recoded")

### Stroke ---------------------------------------------------------------------

temp <- read.table("EUR (main)/Stroke_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz", 
                   header=F, nrows=2)
colnames(temp) <- limma::strsplit2(readLines("EUR (main)/Stroke_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz", n=1),
                                   "\t")
readLines("EUR (main)/Stroke_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz", n=1)
temp

stroke <- data.frame(pheno = "stroke",
                     filename= "Stroke_Bothsex_eur_inv_var_meta_GBMI_052021_nbbkgt1.txt.gz",
                     sep ="\t",
                     snp_col = "rsid",
                     beta_col = "inv_var_meta_beta",
                     se_col = "inv_var_meta_sebeta",
                     effect_allele_col="REF",
                     other_allele_col = "ALT",
                     eaf_col = "all_meta_AF",
                     pval_col = "inv_var_meta_p",
                     note = "beta")


### Type 2 diabetes ------------------------------------------------------------ 

temp <- read.table("EUR (main)/Xue_et_al_T2D_META_Nat_Commun_2018.gz", 
                   header=T, nrows=2)
readLines("EUR (main)/Xue_et_al_T2D_META_Nat_Commun_2018.gz",n=2)
temp

T2D <- data.frame(pheno = "T2D",
                  filename= "Xue_et_al_T2D_META_Nat_Commun_2018.gz",
                  sep =" ",
                  snp_col = "SNP",
                  beta_col = "b",
                  se_col = "se",
                  effect_allele_col="A1",
                  other_allele_col = "A2",
                  eaf_col = "frq_A1",
                  pval_col = "P",
                  note = "beta")

#-------------------------------------------------------------------------------
### UNHEALTHY BEHAVIORS 
#-------------------------------------------------------------------------------

### Alcohol use1 ---------------------------------------------------------------

# consists of two different outcomes, 1) alcohol consumption (C) and 2) alcohol problems (P).

temp <- read.table("EUR (main)/AUDIT_UKB_2018_AJP.txt.gz", 
                   header=T, nrows=2)
temp
readLines("EUR (main)/AUDIT_UKB_2018_AJP.txt.gz",n=1)

auditC <- data.frame(pheno = "auditC",
                     filename= "AUDIT_UKB_2018_AJP.txt.gz",
                     sep =" ",
                     snp_col = "rsid",
                     beta_col = "beta_C",
                     se_col = "se_C",
                     effect_allele_col="a_0",
                     other_allele_col = "a_1",
                     eaf_col = "info",
                     pval_col = "p_C",
                     note = "beta")

auditP <- data.frame(pheno = "auditP",
                     filename= "AUDIT_UKB_2018_AJP.txt.gz",
                     sep =" ",
                     snp_col = "rsid",
                     beta_col = "beta_P",
                     se_col = "se_P",
                     effect_allele_col="a_0",
                     other_allele_col = "a_1",
                     eaf_col = "info",
                     pval_col = "p_P",
                     note = "beta")

### Physical inactivity --------------------------------------------------------

temp <- read.table("EUR (main)/GCST90104341_buildGRCh37.tsv", 
                   header=T, nrows=2)
temp
readLines("EUR (main)/GCST90104341_buildGRCh37.tsv",n=1)

exercise <- data.frame(pheno = "exercise",
                       filename= "GCST90104341_buildGRCh37.tsv",
                       sep ="\t",
                       snp_col = "rsid",
                       beta_col = "BETA",
                       se_col = "SE",
                       effect_allele_col="Allele1",
                       other_allele_col = "Allele2",
                       eaf_col = "Freq1",
                       pval_col = "p_value",
                       note = "beta")

### Smoking --------------------------------------------------------------------

# to recode from VCF format 
temp <- read.table("EUR (main)/Smoking_initiation.vcf.gz", 
                   header=F, nrows=2)
colnames(temp) <- limma::strsplit2(readLines("EUR (main)/Smoking_initiation.vcf.gz", n=110)[110],
                                   "\t")
temp
# smoking <- data.table::fread("EUR (main)/Smoking_initiation.vcf.gz", skip =109) 
# head(smoking[,])
# smoking2 <- cbind(smoking[,1:8],
#                  limma::strsplit2(smoking$`UKB-b-15541`,":"))
# colnames(smoking2)[9:13] <- c(limma::strsplit2(smoking$FORMAT[1],":"))
# smoking2$pval <- 10^-as.numeric(smoking2$LP)
# smoking2 <- smoking2[,-13]
# write.table(smoking2, file = "EUR (main)/Smoking_initiation_recoded.txt",sep="\t",
#             col.names=T, row.names = F, quote=F)
# rm(smoking, smoking2)

smoking <- data.frame(pheno = "smoking",
                      filename= "Smoking_initiation_recoded.txt",
                      sep ="\t",
                      snp_col = "ID",
                      beta_col = "ES",
                      se_col = "SE",
                      effect_allele_col="REF",
                      other_allele_col = "ALT",
                      eaf_col = "AF",
                      pval_col = "pval",
                      note = "recoded")

#-------------------------------------------------------------------------------
### READ OUT
#-------------------------------------------------------------------------------

psych <- rbind(anorexia, asd, bipolar,cannabis,depression2023,anxiety, ADHD2022,
               ocd,ptsd, scz, tourette, suicide2020, opioidDepVexp,
               asthma,CKD, COPD,CAD2022, obesity,stroke,T2D,
               auditC,auditP, exercise, smoking)

colnames(psych)[1] <- "psych"
for(i in 1:ncol(psych)){psych[,i] <- as.character(psych[,i])}

save(psych, file = "EUR (main)/main_outcomes_2023-04-19.Rdata")
load("EUR (main)/main_outcomes_2023-04-19.Rdata")
psych

