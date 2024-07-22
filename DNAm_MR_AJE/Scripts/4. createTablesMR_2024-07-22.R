# A script to make all tables and figures
# Isabel Schuurmans
# 2023-04-20

# please note I inserted scientific p values manually!

# open library
library(ggplot2)
library(cowplot)
library(xlsx)
library(data.table)
library(readr)
library(ggnewscale)

# load in the most recent results
load("YOURPATH/Analysis/Results/compiledMR_age7_2023-05-12.Rdata")
load("YOURPATH/Analysis/Results/compiledMR_age15_2023-05-12.Rdata")
load("YOURPATH/Analysis/Results/compiledMR_marzi_2023-05-12.Rdata")
load("YOURPATH/Analysis/Results/compiledMR_PACE_2023-05-13.Rdata")

# remove outcomes we dont want to focus on
age7results <- age7results[!age7results$id.outcome %in% c("auditT","adhd","depressionPGC_UKB", "suicideAttemptsCorr", "depressionGian","opioid","cardiogram","suicideAttempts","depressionHoward","aud"),] 
age15results <- age15results[!age15results$id.outcome %in% c("auditT","adhd","depressionPGC_UKB", "suicideAttemptsCorr", "depressionGian","opioid", "cardiogram","suicideAttempts","depressionHoward","aud"),]
marziresults <- marziresults[!marziresults$id.outcome %in% c("auditT","adhd","depressionPGC_UKB", "suicideAttemptsCorr", "depressionGian","opioid", "cardiogram","suicideAttempts","depressionHoward","aud"),]
PACEresults <- PACEresults[!PACEresults$id.outcome %in% c("auditT","adhd","depressionPGC_UKB", "suicideAttemptsCorr", "depressionGian","opioid", "cardiogram","suicideAttempts","depressionHoward","aud"),]

# convert physical activity; in order for the outcome have higher values correspond with worse outcomes
age7results[which(age7results$id.outcome == 'exercise'),c("b", "b_adj")] <- age7results[which(age7results$id.outcome == 'exercise'),c("b", "b_adj")] *-1
age15results[which(age15results$id.outcome == 'exercise'),c("b", "b_adj")] <- age15results[which(age15results$id.outcome == 'exercise'),c("b", "b_adj")] *-1
marziresults[which(marziresults$id.outcome == 'exercise'),c("b", "b_adj")] <- marziresults[which(marziresults$id.outcome == 'exercise'),c("b", "b_adj")] *-1
PACEresults[which(PACEresults$id.outcome == 'exercise'),c("b", "b_adj")] <- PACEresults[which(PACEresults$id.outcome == 'exercise'),c("b", "b_adj")] *-1

# load in CPGs
load("C:/Users/isabe/Dropbox (Partners HealthCare)/ALSPAC/Epigenetic/Paper Draft - MR/Analysis/Data/cpgsMR_2023-04-13.Rdata")

# get LL & UL beta
age7results$LL <- age7results$b - 1.96*age7results$se
age7results$UL <- age7results$b + 1.96*age7results$se
age7results$CI <- paste0("[", sprintf("%.2f", round(age7results$LL,2)),", ", sprintf("%.2f", round(age7results$UL,2)), "]")
age7results$CI_beta <- paste0(sprintf("%.2f", round(age7results$b,2)), " ", age7results$CI)
age15results$LL <- age15results$b - 1.96*age15results$se
age15results$UL <- age15results$b + 1.96*age15results$se
age15results$CI <- paste0("[", sprintf("%.2f", round(age15results$LL,2)),", ", sprintf("%.2f", round(age15results$UL,2)), "]")
age15results$CI_beta <- paste0(sprintf("%.2f", round(age15results$b,2)), " ", age15results$CI)
marziresults$LL <- marziresults$b - 1.96*marziresults$se
marziresults$UL <- marziresults$b + 1.96*marziresults$se
marziresults$CI <- paste0("[", sprintf("%.2f", round(marziresults$LL,2)),", ", sprintf("%.2f", round(marziresults$UL,2)), "]")
marziresults$CI_beta <- paste0(sprintf("%.2f", round(marziresults$b,2)), " ", marziresults$CI)
PACEresults$LL <- PACEresults$b - 1.96*PACEresults$se
PACEresults$UL <- PACEresults$b + 1.96*PACEresults$se
PACEresults$CI <- paste0("[", sprintf("%.2f", round(PACEresults$LL,2)),", ", sprintf("%.2f", round(PACEresults$UL,2)), "]")
PACEresults$CI_beta <- paste0(sprintf("%.2f", round(PACEresults$b,2)), " ", PACEresults$CI)

# get odds ratios
age7results$odds[!age7results$outcome %in% c("auditC", "auditP", "auditT")] <- exp(age7results$b)[!age7results$outcome %in% c("auditC", "auditP", "auditT")]
age15results$odds[!age15results$outcome %in% c("auditC", "auditP", "auditT")] <- exp(age15results$b)[!age15results$outcome %in% c("auditC", "auditP", "auditT")]
marziresults$odds[!marziresults$outcome %in% c("auditC", "auditP", "auditT")] <- exp(marziresults$b)[!marziresults$outcome %in% c("auditC", "auditP", "auditT")]
PACEresults$odds[!PACEresults$outcome %in% c("auditC", "auditP", "auditT")] <- exp(PACEresults$b)[!PACEresults$outcome %in% c("auditC", "auditP", "auditT")]

# get CI odds ratios
age7results$CI_odds <- paste0(sprintf("%.2f", round(age7results$odds,2)), " [", 
                              sprintf("%.2f", round(exp(age7results$LL),2)),", ", 
                              sprintf("%.2f", round(exp(age7results$UL),2)), "]")
age15results$CI_odds <- paste0(sprintf("%.2f", round(age15results$odds,2)), " [", 
                              sprintf("%.2f", round(exp(age15results$LL),2)),", ", 
                              sprintf("%.2f", round(exp(age15results$UL),2)), "]")
marziresults$CI_odds <- paste0(sprintf("%.2f", round(marziresults$odds,2)), " [", 
                              sprintf("%.2f", round(exp(marziresults$LL),2)),", ", 
                              sprintf("%.2f", round(exp(marziresults$UL),2)), "]")
PACEresults$CI_odds <- paste0(sprintf("%.2f", round(PACEresults$odds,2)), " [", 
                               sprintf("%.2f", round(exp(PACEresults$LL),2)),", ", 
                               sprintf("%.2f", round(exp(PACEresults$UL),2)), "]")

age7results$CI_odds[is.na(age7results$odds)] <- NA 
age15results$CI_odds[is.na(age15results$odds)] <- NA 
marziresults$CI_odds[is.na(marziresults$odds)] <- NA 
PACEresults$CI_odds[is.na(PACEresults$odds)] <- NA 

# get role DNAm
age7results$roleDNAm <- ifelse(age7results$b_adj > 0, "Risk increasing", "Risk suppressing")
age15results$roleDNAm <- ifelse(age15results$b_adj > 0, "Risk increasing", "Risk suppressing")
marziresults$roleDNAm <- ifelse(marziresults$b_adj > 0, "Risk increasing", "Risk suppressing")
PACEresults$roleDNAm <- ifelse(PACEresults$b_adj > 0, "Risk increasing", "Risk suppressing")

# define categories
psych <- c("anorexia", "asd", "bipolar","cannabis","depression2023","anxiety", "ADHD2022",
  "ocd","ptsd", "scz", "tourette", "suicide2020", "opioidDepVexp")
phys <- c("asthma","CKD", "COPD","CAD2022", "obesity","stroke","T2D")
behav <- c("auditC","auditP", "exercise", "smoking")
  
# select categories
age7results$category <- NA
age7results$category[age7results$outcome %in% psych] <- "Mental health outcomes"
age7results$category[age7results$outcome %in% phys] <- "Physical health outcomes"
age7results$category[age7results$outcome %in% behav] <- "Unhealthy behaviors"

# change order levels for the plot
age7results$category <- factor(age7results$category, levels = c("Mental health outcomes","Physical health outcomes","Unhealthy behaviors"))

# select categories
age15results$category <- NA
age15results$category[age15results$outcome %in% psych] <- "Mental health outcomes"
age15results$category[age15results$outcome %in% phys] <- "Physical health outcomes"
age15results$category[age15results$outcome %in% behav] <- "Unhealthy behaviors"

# change order levels for the plot
age15results$category <- factor(age15results$category, levels = c("Mental health outcomes","Physical health outcomes","Unhealthy behaviors"))

# select categories
marziresults$category <- NA
marziresults$category[marziresults$outcome %in% psych] <- "Mental health outcomes"
marziresults$category[marziresults$outcome %in% phys] <- "Physical health outcomes"
marziresults$category[marziresults$outcome %in% behav] <- "Unhealthy behaviors"

# change order levels for the plot
marziresults$category <- factor(marziresults$category, levels = c("Mental health outcomes","Physical health outcomes","Unhealthy behaviors"))

# select categories
PACEresults$category <- NA
PACEresults$category[PACEresults$outcome %in% psych] <- "Mental health outcomes"
PACEresults$category[PACEresults$outcome %in% phys] <- "Physical health outcomes"
PACEresults$category[PACEresults$outcome %in% behav] <- "Unhealthy behaviors"

# change order levels for the plot
PACEresults$category <- factor(PACEresults$category, levels = c("Mental health outcomes","Physical health outcomes","Unhealthy behaviors"))


#--------------------------------------------------------------------------------------------------------------
# FIGURE S2

# define age
age7results$age <- 'Childhood DNAm'
age15results$age <- 'Adolescent DNAm'
marziresults$age <- "Young adulthood DNAm (triangulation with Marzi et al. (2019))"
PACEresults$age <- "Birth DNAm (triangulation with Ruehlmann et al. (2023))"

# merge
combined_results <- rbind(age7results, age15results)
combined_results <- rbind(combined_results, marziresults)
combined_results <- rbind(combined_results, PACEresults)
combined_results$age <- factor(combined_results$age, levels = c("Childhood DNAm","Adolescent DNAm", 
                                                                "Birth DNAm (triangulation with Ruehlmann et al. (2023))", 
                                                                "Young adulthood DNAm (triangulation with Marzi et al. (2019))"))

# make signficance indicator
combined_results
combined_results$p_sig_fdr <- ifelse(combined_results$fdr < .05, TRUE, NA)
combined_results$p_sig <- ifelse(combined_results$pval < .01 & combined_results$fdr > .05, TRUE, NA)

# select only the result needed anyways
combined_results_fig1 <- combined_results[!is.na(combined_results$fdr) & 
                                            combined_results$id.exposure %in% c(unique(age7results$id.exposure), unique(age15results$id.exposure)),]

# get heatmap order
wide <- reshape(combined_results_fig1[,c("exposure",'outcome',"b_adj")], idvar = "exposure", timevar = "outcome", direction = "wide")
rownames(wide) <- wide[,1]

# psych
wide_psych <- wide[,c(paste0("b_adj.", psych, sep = ""))]
wide_psych <- wide_psych[!rowSums(is.na(wide_psych))>=ncol(wide_psych)-1,]
hm <- heatmap(as.matrix(wide_psych))
psych_order <- names(wide_psych)[hm$colInd]

# phys health
wide_phys <- wide[,c(paste0("b_adj.", phys, sep = ""))]
wide_phys <- wide_phys[!rowSums(is.na(wide_phys))>=ncol(wide_phys)-1,]
hm <- heatmap(as.matrix(wide_phys))
phys_order <- names(wide_phys)[hm$colInd]

# behav
wide_behav <- wide[,c(paste0("b_adj.", behav, sep = ""))]
wide_behav <- wide_behav[!rowSums(is.na(wide_behav))>=ncol(wide_behav)-1,]
hm <- heatmap(as.matrix(wide_behav))
behav_order <- names(wide_behav)[hm$colInd]

# combine names 
order_names <- sub("......", '', c(psych_order, phys_order, behav_order))

# reorder plot outcomes
combined_results_fig1$outcome_plot <- factor(combined_results_fig1$outcome, levels = order_names)

## now also for predictors

# age7
wide_7 <- reshape(combined_results_fig1[combined_results_fig1$age == "Childhood DNAm",c("exposure",'outcome',"b_adj")], idvar = "exposure", timevar = "outcome", direction = "wide")
wide_7 <- wide_7[!rowSums(is.na(wide_7))>=ncol(wide_7)-2,]
hm <- heatmap(as.matrix(wide_7[,-1]))
cpg7_order <- as.character(wide_7[,1][hm$rowInd])
# add missing cpgs
cpg7_order <- c(cpg7_order, unique(age7results$id.exposure)[!unique(age7results$id.exposure) %in% cpg7_order])

# age15
wide_15 <- reshape(combined_results_fig1[combined_results_fig1$age == "Adolescent DNAm",c("exposure",'outcome',"b_adj")], idvar = "exposure", timevar = "outcome", direction = "wide")
hm <- heatmap(as.matrix(wide_15[,-1]))
cpg15_order <- as.character(wide_15[,1][hm$rowInd])
# add missing cpgs
cpg15_order <- c(cpg15_order, unique(age15results$id.exposure)[!unique(age15results$id.exposure) %in% cpg15_order])

# combine names 
order_names <- c(cpg7_order,cpg15_order)

# reorder plot outcomes
combined_results_fig1$exposure_plot <- factor(combined_results_fig1$exposure, levels = order_names)

# rename outcomes in plot -> ALWAYS CHECK THIS PART TO SEE IF NAMES ARE CORRESPONDING!
levels(combined_results_fig1$outcome_plot)

names_outcomes <- c("Depression","Cannabis use disorder", "ADHD",
                    "Suicide attempt",  "Schizophrenia","Bipolar disorder",
                    "Autism spectrum disorder", "PTSD","Opioid dependence",
                    "Anorexia nervosa","Tourette syndrome","Anxiety disorder",
                    "Obsessive-compulsive disorder", 
                    
                    "Stroke", "Chronic kidney disease","Coronary artery disease", 
                    "Type 2 Diabetes" ,"Obesity","Asthma","COPD",
                    
                    "Alcohol problems","Smoking", 
                    "Alcohol consumption",  "Physical inactivity")

combined_results_fig1$outcome_plot_new_names <- combined_results_fig1$outcome_plot
levels(combined_results_fig1$outcome_plot_new_names) <- names_outcomes

levels(combined_results_fig1$category) <- c("Mental health outcomes", "Physical health outcomes", "Unhealthy \nbehaviors")

# plot
plot7 <- ggplot(data = combined_results_fig1, aes(x = exposure_plot, y = outcome_plot_new_names)) +
  geom_tile(aes(fill = b_adj)) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "salmon", limits = c(-1,1), breaks = c(-0.9,0,0.9),
                       labels = c("-1.0  Risk suppressing effects",
                                  " 0.0  ",
                                  " 1.0  Risk increasing effects")) + 
  labs(fill="Effect estimate") +
  theme(
    text = element_text(size = 18),
    panel.background = element_rect(fill = "gray 50"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank())+
  ylab('')+
  xlab('') + 
  theme(axis.text.x = element_text(angle =45, hjust=1)) +
  new_scale("size")+
  geom_point(aes(x = exposure, y = outcome_plot_new_names, size= as.numeric(p_sig)), shape = 20) +
  scale_size(labels = "") + 
  labs(size = "p < .01 & q > .05") +
  new_scale("size")+
  geom_point(aes(x = exposure, y = outcome_plot_new_names, size = as.numeric(p_sig_fdr)), shape = 8) +
  scale_size(labels = c("")) + 
  labs(size = "q < .05") +
  guides(size = guide_legend(order=2))+
  facet_grid(category ~ age, scales = "free", space = 'free')
plot7  

# write out
tiff("YOURPATH/Analysis/FiguresTables/FigureS2-2023-05-04.tiff", width = 7000, height = 4000, res = 400, pointsize = 4, type = "cairo")
plot7
dev.off()

## ----------------------------------------------------------------------------------------

## FIGURE S3

# select only the result needed anyways
combined_results_figS2 <- combined_results[!is.na(combined_results$fdr) & 
                                            combined_results$id.exposure %in% c(unique(marziresults$id.exposure), unique(PACEresults$id.exposure)),]

# get heatmap order
wide <- reshape(combined_results_figS2[,c("exposure",'outcome',"b_adj")], idvar = "exposure", timevar = "outcome", direction = "wide")
rownames(wide) <- wide[,1]

# psych
wide_psych <- wide[,c(paste0("b_adj.", psych, sep = ""))]
wide_psych <- wide_psych[!rowSums(is.na(wide_psych))>=ncol(wide_psych)-1,]
hm <- heatmap(as.matrix(wide_psych))
psych_order <- names(wide_psych)[hm$colInd]

# phys health
wide_phys <- wide[,c(paste0("b_adj.", phys, sep = ""))]
wide_phys <- wide_phys[!rowSums(is.na(wide_phys))>=ncol(wide_phys)-1,]
hm <- heatmap(as.matrix(wide_phys))
phys_order <- names(wide_phys)[hm$colInd]

# behav
wide_behav <- wide[,c(paste0("b_adj.", behav, sep = ""))]
wide_behav <- wide_behav[!rowSums(is.na(wide_behav))>=ncol(wide_behav)-1,]
hm <- heatmap(as.matrix(wide_behav))
behav_order <- names(wide_behav)[hm$colInd]

# combine names 
order_names <- sub("......", '', c(psych_order, phys_order, behav_order))

# reorder plot outcomes
combined_results_figS2$outcome_plot <- factor(combined_results_figS2$outcome, levels = order_names)

## now also for predictors

# PACE
wide_PACE <- reshape(combined_results_figS2[combined_results_figS2$age == "Birth DNAm (triangulation with Ruehlmann et al. (2023))",c("exposure",'outcome',"b_adj")], idvar = "exposure", timevar = "outcome", direction = "wide")
wide_PACE <- wide_PACE[!rowSums(is.na(wide_PACE))>=ncol(wide_PACE)-2,]
hm <- heatmap(as.matrix(wide_PACE[,-1]))
cpgPACE_order <- as.character(wide_PACE[,1][hm$rowInd])
# add missing cpgs
cpgPACE_order <- c(cpgPACE_order, unique(PACEresults$id.exposure)[!unique(PACEresults$id.exposure) %in% cpgPACE_order])

# marzi
wide_marzi <- reshape(combined_results_figS2[combined_results_figS2$age == "Young adulthood DNAm (triangulation with Marzi et al. (2019))",c("exposure",'outcome',"b_adj")], idvar = "exposure", timevar = "outcome", direction = "wide")
hm <- heatmap(as.matrix(wide_marzi[,-1]))
cpgmarzi_order <- as.character(wide_marzi[,1][hm$rowInd])
# add missing cpgs
cpgmarzi_order <- c(cpgmarzi_order, unique(marziresults$id.exposure)[!unique(marziresults$id.exposure) %in% cpgmarzi_order])

# combine names 
order_names <- c(cpgPACE_order,cpgmarzi_order)

# reorder plot outcomes
combined_results_figS2$exposure_plot <- factor(combined_results_figS2$exposure, levels = order_names)

# rename outcomes in plot -> ALWAYS CHECK THIS PART TO SEE IF NAMES ARE CORRESPONDING!
levels(combined_results_figS2$outcome_plot)

names_outcomes <- c("Opioid dependence","Bipolar disorder","Schizophrenia",
                    "Suicide attempt", "Depression","Cannabis use disorder", "ADHD",
                    "Autism spectrum disorder", "PTSD",
                    "Anorexia nervosa","Tourette syndrome","Anxiety disorder",
                    "Obsessive-compulsive disorder", 
                    
                    "Type 2 Diabetes" ,"Coronary artery disease", 
                    "Stroke", "COPD","Chronic kidney disease","Obesity","Asthma",
                    
                    
                    "Alcohol problems","Alcohol consumption", "Smoking", 
                     "Physical inactivity")

combined_results_figS2$outcome_plot_new_names <- combined_results_figS2$outcome_plot
levels(combined_results_figS2$outcome_plot_new_names) <- names_outcomes

# fix labels
swr = function(string, nwrap=30) {
  paste(strwrap(string, width=nwrap), collapse="\n")
}
swr = Vectorize(swr)

# Create line breaks in Year
combined_results_figS2$age_n = as.factor(swr(combined_results_figS2$age))
combined_results_figS2$age_n <- factor(combined_results_figS2$age_n, levels = rev(levels(as.factor(combined_results_figS2$age_n))))

levels(combined_results_figS2$category) <- c("Mental health outcomes", "Physical health outcomes", "Unhealthy \nbehaviors")

# plot
plot7 <- ggplot(data = combined_results_figS2, aes(x = exposure_plot, y = outcome_plot_new_names)) +
  geom_tile(aes(fill = b_adj)) +
  scale_fill_gradient2(low = "steelblue", mid = "white", high = "salmon", limits = c(-1,1), breaks = c(-0.9,0,0.9),
                       labels = c("-1.0  Risk suppressing effects",
                                  " 0.0  ",
                                  " 1.0  Risk increasing effects")) + 
  labs(fill="Effect estimate") +
  theme(
    text = element_text(size = 18),
    panel.background = element_rect(fill = "gray 50"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank())+
  ylab('')+
  xlab('') + 
  theme(axis.text.x = element_text(angle =45, hjust=1)) +
  new_scale("size")+
  geom_point(aes(x = exposure, y = outcome_plot_new_names, size= as.numeric(p_sig)), shape = 20) +
  scale_size(labels = "") + 
  labs(size = "p < .01 & q > .05") +
  new_scale("size")+
  geom_point(aes(x = exposure, y = outcome_plot_new_names, size = as.numeric(p_sig_fdr)), shape = 8) +
  scale_size(labels = c("")) + 
  labs(size = "q < .05") +
  guides(size = guide_legend(order=2))+
  facet_grid(category ~ age_n, scales = "free", space = 'free')
plot7  

# write out
tiff("YOURPATH/Analysis/FiguresTables/FigureS3-2023-05-04.tiff", width = 7000, height = 4000, res = 400, pointsize = 4, type = "cairo")
plot7
dev.off()

## ----------------------------------------------------------------------------------------
## TABLE 2

# Top hits age 7

# get ones that are significant
age7results_plot <- age7results[which(age7results$pval < .01 & !is.na(age7results$fdr)),]

# select columns
tabl2_age7 <- age7results_plot[,c("category","outcome", "exposure", "method","nsnp","CI_beta",'CI_odds', "pval", "fdr", "roleDNAm")]

# round up
tabl2_age7$pval[tabl2_age7$pval > 0.0005] <-  as.character(paste0(sprintf("%.3f", round(tabl2_age7$pval[tabl2_age7$pval > 0.0005],3))))
tabl2_age7$fdr[tabl2_age7$fdr > 0.0005] <-  as.character(paste0(sprintf("%.3f", round(tabl2_age7$fdr[tabl2_age7$fdr > 0.0005],3))))

ifelse(age7results_plot$b/age7results_plot$b_adj == -1, "negative", "positive")

# write out
write.xlsx(tabl2_age7, "YOURPATH/Analysis/FiguresTables/Table2-2023-05-04.xlsx")

## ----------------------------------------------------------------------------------------
## TABLE 3

# Top hits age 15

# get ones that are significant
age15results_plot <- age15results[which(age15results$pval < .01 & !is.na(age15results$fdr)),]

# select columns
tabl2_age15 <- age15results_plot[,c("category","outcome", "exposure", "method","nsnp","CI_beta",'CI_odds', "pval", "fdr", "roleDNAm")]

# round up
tabl2_age15$pval[tabl2_age15$pval > 0.0005] <-  as.character(paste0(sprintf("%.3f", round(tabl2_age15$pval[tabl2_age15$pval > 0.0005],3))))
tabl2_age15$fdr[tabl2_age15$fdr > 0.0005]  <-  as.character(paste0(sprintf("%.3f", round(tabl2_age15$fdr[tabl2_age15$fdr > 0.0005],3))))

write.xlsx(tabl2_age15, "YOURPATH/Analysis/FiguresTables/Table3-2023-05-04.xlsx")

## ----------------------------------------------------------------------------------------
## FIGURE 1

# define ages data
tabl2_age7$age <- "Childhood DNAm (primary analyses based on Lussier et al. (2019))"
tabl2_age15$age <- "Adolescent DNAm (primary analyses based on Lussier et al. (2023))"
tabl2_repli_birth <- tabl2_repli[tabl2_repli$age == "Birth DNAm (triangulation with Ruehlmann et al. (2023))",]
table2_repli_adol <- tabl2_repli[tabl2_repli$age == "Young adulthood DNAm (triangulation with Marzi et al. (2019))",]

# bind
#allresults <- rbind(tabl2_age7, tabl2_age15, tabl2_repli)
#allresults

# get number of associations per outcome
age7 <- data.frame(table(tabl2_age7$outcome))
age15 <- data.frame(table(tabl2_age15$outcome))
birth <- data.frame(table(tabl2_repli$outcome[tabl2_repli$age == "Birth DNAm (triangulation with Ruehlmann et al. (2023))"]))
adol <- data.frame(table(tabl2_repli$outcome[tabl2_repli$age == "Young adulthood DNAm (triangulation with Marzi et al. (2019))"]))

# get role dna m age 7
age7$rolednam <- NA
tabl2_age7$roleDNAm_num <- ifelse(tabl2_age7$roleDNAm == "Risk increasing", 1, 2)

for (i in 1:nrow(age7)){
  
  if (age7[i,"Freq"] == 0){
    age7$rolednam[i] <- 0
  }
  
  else {
    
    if (age7[i,"Freq"] == 1){
      age7$rolednam[i] <- tabl2_age7[tabl2_age7$outcome == age7[i,"Var1"], "roleDNAm_num"] 
    } 
    
    else {
      
      age7$rolednam[i] <- ifelse(
        var(tabl2_age7[tabl2_age7$outcome == age7[i,"Var1"], "roleDNAm_num"]) == 0,
        tabl2_age7[tabl2_age7$outcome == age7[i,"Var1"], "roleDNAm_num"][1],3)
    }
  }
}

# get role dna m age 15
age15$rolednam <- NA
tabl2_age15$roleDNAm_num <- ifelse(tabl2_age15$roleDNAm == "Risk increasing", 1, 2)

for (i in 1:nrow(age15)){
  
  if (age15[i,"Freq"] == 0){
    age15$rolednam[i] <- 0
  }
  
  else {
    
    if (age15[i,"Freq"] == 1){
      age15$rolednam[i] <- tabl2_age15[tabl2_age15$outcome == age15[i,"Var1"], "roleDNAm_num"] 
    } 
    
    else {
      
      age15$rolednam[i] <- ifelse(
        var(tabl2_age15[tabl2_age15$outcome == age15[i,"Var1"], "roleDNAm_num"]) == 0,
        tabl2_age15[tabl2_age15$outcome == age15[i,"Var1"], "roleDNAm_num"][1],3)
    }
  }
}

# get role dna m repli birth
birth$rolednam <- NA
tabl2_repli_birth$roleDNAm_num <- ifelse(tabl2_repli_birth$roleDNAm == "Risk increasing", 1, 2)

for (i in 1:nrow(birth)){
  
  if (birth[i,"Freq"] == 0){
    birth$rolednam[i] <- 0
  }
  
  else {
    
    if (birth[i,"Freq"] == 1){
      birth$rolednam[i] <- tabl2_repli_birth[tabl2_repli_birth$outcome == birth[i,"Var1"], "roleDNAm_num"] 
    } 
    
    else {
      
      birth$rolednam[i] <- ifelse(
        var(tabl2_repli_birth[tabl2_repli_birth$outcome == birth[i,"Var1"], "roleDNAm_num"]) == 0,
        tabl2_repli_birth[tabl2_repli_birth$outcome == birth[i,"Var1"], "roleDNAm_num"][1],3)
    }
  }
}

# get role dna m repli adol
adol$rolednam <- NA
table2_repli_adol$roleDNAm_num <- ifelse(table2_repli_adol$roleDNAm == "Risk increasing", 1, 2)

for (i in 1:nrow(adol)){
  
  if (adol[i,"Freq"] == 0){
    adol$rolednam[i] <- 0
  }
  
  else {
    
    if (adol[i,"Freq"] == 1){
      adol$rolednam[i] <- table2_repli_adol[table2_repli_adol$outcome == adol[i,"Var1"], "roleDNAm_num"] 
    } 
    
    else {
      
      adol$rolednam[i] <- ifelse(
        var(table2_repli_adol[table2_repli_adol$outcome == adol[i,"Var1"], "roleDNAm_num"]) == 0,
        table2_repli_adol[table2_repli_adol$outcome == adol[i,"Var1"], "roleDNAm_num"][1],3)
    }
  }
}


# bind number of associations 
newplot <- data.frame(unique(age7results$outcome))
newplot <- merge(newplot, birth, by.x = "unique.age7results.outcome.", by.y = 'Var1')
newplot <- merge(newplot, age7, by.x = "unique.age7results.outcome.", by.y = 'Var1')
newplot <- merge(newplot, age15, by.x = "unique.age7results.outcome.", by.y = 'Var1')
newplot <- merge(newplot, adol, by.x = "unique.age7results.outcome.", by.y = 'Var1')

# name them      
names(newplot) <- c("outcome", "birth", "age7", "age15", "age18")

# now melt the data
melt_plot <- data.frame(rep(c("birth","age7","age15","age18"), each = nrow(newplot)))
melt_plot$outcomes <- rep(newplot$outcome, times = 4)
melt_plot$number_asso <- c(newplot[,2], newplot[,4], newplot[,6], newplot[,8])
melt_plot$direction <- c(newplot[,3], newplot[,5], newplot[,7], newplot[,9])
names(melt_plot) <- c("time","outcomes", "times","direction")

# names direction
melt_plot$direction <- as.factor(melt_plot$direction)
levels(melt_plot$direction) <- c("","Increasing risk","Suppressing risk", "Mixed")

# also order outcome and time
melt_plot$time[melt_plot$time == "birth"] <- "Birth (8 DNAm loci)"
melt_plot$time[melt_plot$time == "age7"] <- "Childhood (21 DNAm loci)"
melt_plot$time[melt_plot$time == "age15"] <- "Adolescence (15 DNAm loci)"
melt_plot$time[melt_plot$time == "age18"] <- "Young adulthood (19 DNAm loci)"

melt_plot$time <- as.factor(melt_plot$time)
melt_plot$time <- factor(melt_plot$time, levels = c("Birth (8 DNAm loci)", "Childhood (21 DNAm loci)", 
                                                    "Adolescence (15 DNAm loci)", "Young adulthood (19 DNAm loci)"))

# get category
melt_plot$category <- NA

for (i in 1:nrow(melt_plot)) {
  melt_plot$category[i] <- as.character(age7results[age7results$outcome == melt_plot[i,"outcomes"], "category"])[1] }

melt_plot$category[melt_plot$category == "Mental health outcomes"] <- "Mental health"
melt_plot$category[melt_plot$category == "Physical health outcomes"] <- "Physical health"

# rename category
swr = function(string, nwrap=15) {
  paste(strwrap(string, width=nwrap), collapse="\n")
}
swr = Vectorize(swr)

# Create line breaks in Year
melt_plot$category = as.factor(swr(melt_plot$category))


# fix names
melt_plot$outcomes <- as.character(melt_plot$outcomes)
melt_plot[melt_plot$outcomes == "ADHD2022","outcomes"] <- "ADHD"
melt_plot[melt_plot$outcomes == "anorexia","outcomes"] <- "Anorexia nervosa"
melt_plot[melt_plot$outcomes == "asd","outcomes"] <- "Autism spectrum disorder"
melt_plot[melt_plot$outcomes == "asthma","outcomes"] <- "Asthma"
melt_plot[melt_plot$outcomes == "auditC","outcomes"] <- "Alcohol consumption"
melt_plot[melt_plot$outcomes == "auditP","outcomes"] <- "Alcohol problems"
melt_plot[melt_plot$outcomes == "bipolar","outcomes"] <- "Bipolar disorder"
melt_plot[melt_plot$outcomes == "cannabis","outcomes"] <- "Cannabis use disorder"
melt_plot[melt_plot$outcomes == "CKD","outcomes"] <- "Chronic kidney disease"
melt_plot[melt_plot$outcomes == "depression2023","outcomes"] <- "Depression"
melt_plot[melt_plot$outcomes == "COPD","outcomes"] <- "COPD"
melt_plot[melt_plot$outcomes == "stroke","outcomes"] <- "Stroke"
melt_plot[melt_plot$outcomes == "scz","outcomes"] <- "Schizophrenia"
melt_plot[melt_plot$outcomes == "tourette","outcomes"] <- "Tourette syndrome"
melt_plot[melt_plot$outcomes == "ocd","outcomes"] <- "Obsessive-compulsive disorder"
melt_plot[melt_plot$outcomes == "anxiety","outcomes"] <- "Anxiety disorders"
melt_plot[melt_plot$outcomes == "exercise","outcomes"] <- "Physical inactivity"
melt_plot[melt_plot$outcomes == "CAD2022","outcomes"] <- "Coronary artery disease"
melt_plot[melt_plot$outcomes == "suicide2020","outcomes"] <- "Suicide attempts"
melt_plot[melt_plot$outcomes == "ptsd","outcomes"] <- "PTSD"
melt_plot[melt_plot$outcomes == "opioidDepVexp","outcomes"] <- "Opioid dependence"
melt_plot[melt_plot$outcomes == "T2D","outcomes"] <- "Type 2 diabetes"
melt_plot[melt_plot$outcomes == "smoking","outcomes"] <- "Smoking"
melt_plot[melt_plot$outcomes == "obesity","outcomes"] <- "Obesity"



order_factors <- c("Autism spectrum disorder" , "Cannabis use disorder" ,  "Opioid dependence" , 
                   
                   "Anorexia nervosa","ADHD",  "Depression"  , "Suicide attempts"  ,    
                   
                   "Bipolar disorder",  "Schizophrenia" ,   "Obsessive-compulsive disorder", 
                   "PTSD","Anxiety disorders","Tourette syndrome" ,
                   
                   "Type 2 diabetes" ,"Asthma" , "Stroke" ,"COPD" , "Coronary artery disease", "Chronic kidney disease" , 
                   "Obesity", 
                   
                   "Physical inactivity","Alcohol consumption" ,"Alcohol problems" ,"Smoking"
)


order_factors <- order_factors[rev(order(order_factors))]

melt_plot$outcomes <- factor(melt_plot$outcomes, 
                             levels = order_factors)

# define colors
colors <- c("white","indianred2", "skyblue2", "navajowhite")

melt_plot$numbers <- ifelse(melt_plot$times > 1, melt_plot$times, NA)

# plot
plot7 <- ggplot(data = melt_plot, aes(x = time, y = outcomes, fill = factor(direction))) +
  geom_tile(color = 'black') +
  scale_fill_manual(values=colors) +
  labs(fill="") +
  geom_text(aes(label = numbers)) +
  theme_classic()+
  theme(
    text = element_text(size = 18)) +
  ylab('')+
  xlab('') + 
  theme(axis.text.x = element_text(angle =45, hjust=1)) +
  facet_grid(category ~., scales = "free", space = 'free')
plot7  

# write out
tiff("YOURPATH/Analysis/FiguresTables/Figure1-2023-05-04.tiff", 
     width = 2900, height = 4000, res = 400, pointsize = 4, type = "cairo")
plot7
dev.off()



## ----------------------------------------------------------------------------------------
## TABLE 4

# Top hits replication 

# get ones that are significant
marziresults_plot <- marziresults[which(marziresults$pval < .01 & !is.na(marziresults$fdr)),]
PACEresults_plot <- PACEresults[which(PACEresults$pval < .01 & !is.na(PACEresults$fdr)),]

# bind
repliresults <- rbind(PACEresults_plot, marziresults_plot)

# select columns
tabl2_repli <- repliresults[,c("age","category","outcome", "exposure", "method","nsnp","CI_beta",'CI_odds', "pval", "fdr", "roleDNAm")]

# round up
tabl2_repli$pval[tabl2_repli$pval > 0.0005] <-  as.character(paste0(sprintf("%.3f", round(tabl2_repli$pval[tabl2_repli$pval > 0.0005],3))))
tabl2_repli$fdr[tabl2_repli$fdr > 0.0005] <-  as.character(paste0(sprintf("%.3f", round(tabl2_repli$fdr[tabl2_repli$fdr > 0.0005],3))))

write.xlsx(tabl2_repli, "YOURPATH/Analysis/FiguresTables/Table4-2023-05-04.xlsx")

## ----------------------------------------------------------------------------------------
# NUMBERS FOR RESULTS

## AGE 7

# DNAm loci SNPs
name <- unique(age7results$id.exposure)
cpgs_7 <- data.frame(name)
cpgs_7$nsps_max <- NA
for (i in 1:length(name)){
  cpgs_7$nsps_max[i] <- max(age7results[age7results$id.exposure == name[i],"nsnp"])
}
table(cpgs_7$nsps_max)
sum(length(unique(cpgs_7$name)))

# number of associations
length(tabl2_age7$outcome[as.numeric(tabl2_age7$pval)<=.01])

length(tabl2_age7$outcome[as.numeric(tabl2_age7$pval)<=.01 & tabl2_age7$category == "Mental health outcomes"])
(tabl2_age7$outcome[as.numeric(tabl2_age7$pval)<=.01 & tabl2_age7$category == "Mental health outcomes"])

length(tabl2_age7$outcome[as.numeric(tabl2_age7$pval)<=.01 & tabl2_age7$category == "Physical health outcomes"])
(tabl2_age7$outcome[as.numeric(tabl2_age7$pval)<=.01 & tabl2_age7$category == "Physical health outcomes"])

length(tabl2_age7$outcome[as.numeric(tabl2_age7$pval)<=.01 & tabl2_age7$category == "Unhealthy behaviors"])
(tabl2_age7$outcome[as.numeric(tabl2_age7$pval)<=.01 & tabl2_age7$category == "Unhealthy behaviors"])

# results fdr < .05
tabl2_age7[as.numeric(tabl2_age7$fdr)<=.05,c('exposure', 'outcome', 'CI_odds','fdr')]
nrow(tabl2_age7[as.numeric(tabl2_age7$fdr)<=.05,c('exposure', 'outcome', 'CI_odds','fdr')])

# role DNAm
prop.table(table(tabl2_age7[,c("roleDNAm")]))
tabl2_age7

## AGE 15

# DNAm loci SNPs
name <- unique(age15results$id.exposure)
cpgs_15 <- data.frame(name)
cpgs_15$nsps_max <- NA
for (i in 1:length(name)){
  cpgs_15$nsps_max[i] <- max(age15results[age15results$id.exposure == name[i],"nsnp"])
} 
table(cpgs_15$nsps_max)
sum(length(unique(cpgs_15$name)))

# number of associations
length(tabl2_age15$outcome[as.numeric(tabl2_age15$pval)<=.01])

length(tabl2_age15$outcome[as.numeric(tabl2_age15$pval)<=.01 & tabl2_age15$category == "Mental health outcomes"])
(tabl2_age15$outcome[as.numeric(tabl2_age15$pval)<=.01 & tabl2_age15$category == "Mental health outcomes"])

length(tabl2_age15$outcome[as.numeric(tabl2_age15$pval)<=.01 & tabl2_age15$category == "Physical health outcomes"])
(tabl2_age15$outcome[as.numeric(tabl2_age15$pval)<=.01 & tabl2_age15$category == "Physical health outcomes"])

length(tabl2_age15$outcome[as.numeric(tabl2_age15$pval)<=.01 & tabl2_age15$category == "Unhealthy behaviors"])
(tabl2_age15$outcome[as.numeric(tabl2_age15$pval)<=.01 & tabl2_age15$category == "Unhealthy behaviors"])

# results fdr < .05
tabl2_age15[as.numeric(tabl2_age15$fdr)<=.05,c('exposure', 'outcome', 'CI_odds','fdr')]
nrow(tabl2_age15[as.numeric(tabl2_age15$fdr)<=.05,c('exposure', 'outcome', 'CI_odds','fdr')])

# role DNAm
prop.table(table(tabl2_age15[,c("roleDNAm")]))
tabl2_age15

## PACE

# DNAm loci SNPs
name <- unique(PACEresults$id.exposure)
cpgs_PACE <- data.frame(name)
cpgs_PACE$nsps_max <- NA
for (i in 1:length(name)){
  cpgs_PACE$nsps_max[i] <- max(PACEresults[PACEresults$exposure == name[i],"nsnp"])
} 
table(cpgs_PACE$nsps_max)
sum(length(unique(cpgs_PACE$name)))

# number of associations
length(PACEresults_plot$outcome[PACEresults_plot$pval<=.01])
(PACEresults_plot$outcome[PACEresults_plot$pval<=.01])

# results fdr < .05
PACEresults_plot[PACEresults_plot$fdr<=.05,c('exposure', 'outcome', 'CI_odds','fdr')]
nrow(PACEresults_plot[PACEresults_plot$fdr<=.05,c('exposure', 'outcome', 'CI_odds','fdr')])

# role DNAm
prop.table(table(PACEresults_plot[,c("roleDNAm")]))
PACEresults_plot


## MARZI

# DNAm loci SNPs
name <- unique(marziresults$id.exposure)
cpgs_marzi <- data.frame(name)
cpgs_marzi$nsps_max <- NA
for (i in 1:length(name)){
  cpgs_marzi$nsps_max[i] <- max(marziresults[marziresults$exposure == name[i],"nsnp"])
} 
table(cpgs_marzi$nsps_max)
sum(length(unique(cpgs_marzi$name)))

# number of associations
length(tabl2_marzi$outcome[tabl2_marzi$pval<=.01])
(tabl2_marzi$outcome[tabl2_marzi$pval<=.01])

# results fdr < .05
tabl2_marzi[tabl2_marzi$fdr<=.05,c('exposure', 'outcome', 'CI_odds','fdr')]
nrow(tabl2_marzi[tabl2_marzi$fdr<=.05,c('exposure', 'outcome', 'CI_odds','fdr')])

# role DNAm
prop.table(table(tabl2_marzi[,c("roleDNAm")]))
tabl2_marzi


### CHI SQUARE
df1 <- data.frame(table(tabl2_repli[tabl2_repli$age == "Birth DNAm (triangulation with Ruehlmann et al. (2023))",c("roleDNAm")]),
           table(tabl2_age7[,c("roleDNAm")]),
           table(tabl2_age15[,c("roleDNAm")]),
           table(tabl2_repli[tabl2_repli$age == "Young adulthood DNAm (triangulation with Marzi et al. (2019))",c("roleDNAm")]))


df1 <- df1[,c(2,4,6,8)]

chisq.test(df1)
#---------------- SUPPLEMENTAL INFORMATION ------------------------------------------------

## ----------------------------------------------------------------------------------------
## TABLE S1

# Table S1. Annotated DNAm loci

# select information
tableS1 <- cpgsMR[,c("CpG","Age",  "UCSC_gene", "UCSC_location", "Closest_TSS","Adversity", "beta")]

# remove the mediation snps
tableS1 <- tableS1[tableS1$Age != "Mediation",] 

# round beta
tableS1$beta <- round(tableS1$beta, 2)

# fix names age
tableS1$Age[tableS1$Age == "7"] <- "Childhood DNAm"
tableS1$Age[tableS1$Age == "15"] <- "Adolescent DNAm"

# fix beta
tableS1$beta <-  as.character(paste0(sprintf("%.2f", round(tableS1$beta,2))))

# reorder
tableS1 <- tableS1[,c("CpG","Age", "UCSC_gene", "UCSC_location", "Closest_TSS","Adversity", "beta")]
write.xlsx(tableS1, "YOURPATH/Analysis/FiguresTables/TableS1-2023-04-24.xlsx")

## ----------------------------------------------------------------------------------------
## TABLE S2

# Table S2. Annotated instrumental variables

# read in information
load("YOURPATH/Analysis/Results/age7_loci_MRclump_2023-04-18A/exposureSNPs_goDMC_2023-04-18A.Rdata")
exposureSNPs_7 <- exposureSNPs
exposureSNPs_7$age <- "Childhood DNAm" 
load("YOURPATH/Analysis/Results/age15_loci_MRclump_2023-04-19A/exposureSNPs_goDMC_2023-04-19A.Rdata")
exposureSNPs_15 <- exposureSNPs
exposureSNPs_15$age <- "Adolescent DNAm"

# make it a table
exposureSNPs <- rbind(exposureSNPs_7, exposureSNPs_15)

# select columns
tableS2 <- exposureSNPs[,c("age", "cpg","SNP",'snp', "cistrans")]
write.xlsx(tableS2, "YOURPATH/Analysis/FiguresTables/TableS2-2023-04-24.xlsx")

## ----------------------------------------------------------------------------------------
## TABLE S3

#Table S3. Number of DNAm loci analyzed for each health outcome

cpgs_outcome7 <- cpgs_outcome15 <- outcomes <- unique(age7results$id.outcome)
for (i in 1:length(outcomes)) cpgs_outcome7[i] <- length(unique(age7results[age7results$id.outcome == outcomes[i],'id.exposure']))
for (i in 1:length(outcomes)) cpgs_outcome15[i] <- length(unique(age15results[age15results$id.outcome == outcomes[i],'id.exposure']))

cpgs_n <- data.frame(outcomes, cpgs_outcome7, cpgs_outcome15)
write.xlsx(cpgs_n, "YOURPATH/Analysis/FiguresTables/TableS3-2023-04-24.xlsx")

## ----------------------------------------------------------------------------------------
## TABLE S4

# Table S4. To adversity-related DNAm loci in causal relation to negative health outcomes (full results)

# select per time point what is needed
tableS3_7 <- age7results[,c("age",'category', "outcome", "exposure", "method", "nsnp", "CI_beta","CI_odds", "pval", "fdr", "roleDNAm")]
tableS3_15 <- age15results[,c("age",'category', "outcome", "exposure", "method", "nsnp", "CI_beta","CI_odds","pval", "fdr", "roleDNAm")]

# bind 
tableS3 <- rbind(tableS3_7, tableS3_15)

# clean up
tableS3 <- tableS3[!is.na(tableS3$fdr),]
tableS3$pval[tableS3$pval > 0.0005] <-  as.character(paste0(sprintf("%.3f", round(tableS3$pval[tableS3$pval > 0.0005],3))))
tableS3$fdr[tableS3$fdr > 0.0005] <-  as.character(paste0(sprintf("%.3f", round(tableS3$fdr[tableS3$fdr > 0.0005],3))))

# write out
write.xlsx(tableS3, "YOURPATH/Analysis/FiguresTables/TableS4-2023-04-24.xlsx")


## ----------------------------------------------------------------------------------------
## TABLE S5

#Table S5. To adversity-related adolescent DNAm loci in causal relation to health outcomes (full results replication analyses)

# select columns
tables5_1 <- marziresults[,c("age","category","outcome", "exposure", "method","nsnp","CI_beta",'CI_odds', "pval", "fdr", "roleDNAm")]
tables5_2 <- PACEresults[,c("age","category","outcome", "exposure", "method","nsnp","CI_beta",'CI_odds', "pval", "fdr", "roleDNAm")]

# bind
tables5 <- rbind(tables5_1, tables5_2)

# clean up
tables5 <- tables5[!is.na(tables5$fdr),]
tables5$pval[tables5$pval > 0.0005] <-  as.character(paste0(sprintf("%.3f", round(tables5$pval[tables5$pval > 0.0005],3))))
tables5$fdr[tables5$fdr > 0.0005] <-  as.character(paste0(sprintf("%.3f", round(tables5$fdr[tables5$fdr > 0.0005],3))))

write.xlsx(tables5, "YOURPATH/Analysis/FiguresTables/TableS5-2023-05-04.xlsx")

## ----------------------------------------------------------------------------------------
# TABLE S6

# Table S6. Follow-up of DNAm loci with multiple SNPs, using Wald ratio's on all single SNP

# retrieve significant CpGs
age7single$combinedoutcome <- paste0(age7single$outcome, " (", age7single$id.exposure, ")")
age7results_plot$combinedoutcome <- paste0(age7results_plot$outcome, " (", age7results_plot$id.exposure, ")")
tableS6_7 <- age7single[age7single$combinedoutcome %in% age7results_plot$combinedoutcome,]

age15single$combinedoutcome <- paste0(age15single$outcome, " (", age15single$id.exposure, ")")
age15results_plot$combinedoutcome <- paste0(age15results_plot$outcome, " (", age15results_plot$id.exposure, ")")
tableS6_15 <- age15single[age15single$combinedoutcome %in% age15results_plot$combinedoutcome,]

marzisingle$combinedoutcome <- paste0(marzisingle$outcome, " (", marzisingle$id.exposure, ")")
marziresults_plot$combinedoutcome <- paste0(marziresults_plot$outcome, " (", marziresults_plot$id.exposure, ")")
tableS6_m <- marzisingle[marzisingle$combinedoutcome %in% marziresults_plot$combinedoutcome,]

PACEsingle$combinedoutcome <- paste0(PACEsingle$outcome, " (", PACEsingle$id.exposure, ")")
PACEresults_plot$combinedoutcome <- paste0(PACEresults_plot$outcome, " (", PACEresults_plot$id.exposure, ")")
tableS6_p <- PACEsingle[PACEsingle$combinedoutcome %in% PACEresults_plot$combinedoutcome,]

# bind
tableS6_7$age <- "Childhood DNAm"
tableS6_15$age <- "Adolescence DNAm"
tableS6_m$age <- "Young adulthood DNAm (triangulation with Marzi et al. (2019))"
tableS6_p$age  <- "Birth DNAm (triangulation with Ruehlmann et al. (2023))"

tableS6 <- rbind(tableS6_7, tableS6_15)
tableS6 <- rbind(tableS6, tableS6_m)
tableS6 <- rbind(tableS6, tableS6_p)

# remove this with missing p-values
tableS6 <- tableS6[!is.na(tableS6$p),]
# remove also anything that is not single snps
tableS6 <- tableS6[tableS6$SNP != "All - Inverse variance weighted",]
tableS6 <- tableS6[tableS6$SNP != "All - MR Egger",]

# get CIs
tableS6$LL <- tableS6$b - 1.96*tableS6$se
tableS6$UL <- tableS6$b + 1.96*tableS6$se
tableS6$CI <- paste0("[", sprintf("%.2f", round(tableS6$LL,3)),", ", sprintf("%.2f", round(tableS6$UL,2)), "]")
tableS6$CI_beta <- paste0(sprintf("%.2f", round(tableS6$b,2)), " ", tableS6$CI)

# get odds ratios
tableS6$odds[!tableS6$outcome %in% c("auditC", "auditP", "auditT")] <- 
  exp(tableS6$b)[!tableS6$outcome %in% c("auditC", "auditP", "auditT")]

# get CI odds ratios
tableS6$CI_odds <- paste0(sprintf("%.2f", round(tableS6$odds,2)), " [", 
                          sprintf("%.2f", round(exp(tableS6$LL),2)),", ", 
                          sprintf("%.2f", round(exp(tableS6$UL),2)), "]")
tableS6$CI_odds[is.na(tableS6$odds)] <- NA 

# define categories
tableS6$category <- NA
tableS6$category[tableS6$outcome %in% c("anxiety","ADHD2022","anorexia", "asd","aud", "bipolar","cannabis","depression2023"
                                        ,"depressionGian", "ocd","opioid","ptsd", "scz","suicideAttempt", "tourette", "suicide2020")] <- "Mental health outcomes"
tableS6$category[tableS6$outcome %in% c("asthma","CKD", "COPD","cardiogram", "obesity","stroke","T2D", "CAD2022")] <- "Physical health outcomes"
tableS6$category[tableS6$outcome %in% c("auditC","auditP", "auditT","exercise", "smoking")] <- "Unhealthy behaviors"

tableS6[tableS6$p<.001,]

# select and clean up
tableS6 <- tableS6[,c("age","category","id.outcome", "id.exposure", "SNP", "CI_beta", "CI_odds", "p")]
tableS6$p[tableS6$p > .0005] <-  as.character(paste0(sprintf("%.3f", round(tableS6$p[tableS6$p > .0005],3))))

# write out
write.xlsx(tableS6, "YOURPATH/Analysis/FiguresTables/TableS6-2023-04-24.xlsx")

## ----------------------------------------------------------------------------------------
# TABLE S7

# Table S7. Follow-up of DNAm loci with multiple SNPs, using Wald ratio's on all single SNP

# get unique indicator for cpg x outcome
age7pleio$combinedoutcome <- paste0(age7pleio$outcome, " (", age7pleio$id.exposure, ")")
age7results_plot$combinedoutcome <- paste0(age7results_plot$outcome, " (", age7results_plot$id.exposure, ")")
age15pleio$combinedoutcome <- paste0(age15pleio$outcome, " (", age15pleio$id.exposure, ")")
age15results_plot$combinedoutcome <- paste0(age15results_plot$outcome, " (", age15results_plot$id.exposure, ")")
marzipleio$combinedoutcome <- paste0(marzipleio$outcome, " (", marzipleio$id.exposure, ")")
marziresults_plot$combinedoutcome <- paste0(marziresults_plot$outcome, " (", marziresults_plot$id.exposure, ")")
PACEpleio$combinedoutcome <- paste0(PACEpleio$outcome, " (", PACEpleio$id.exposure, ")")
PACEresults_plot$combinedoutcome <- paste0(PACEresults_plot$outcome, " (", PACEresults_plot$id.exposure, ")")

# retrieve significant CpGs
tableS4_7 <- age7pleio[age7pleio$combinedoutcome %in% age7results_plot$combinedoutcome,]
tableS4_15 <- age15pleio[age15pleio$combinedoutcome %in% age15results_plot$combinedoutcome,]
tableS4_m <- marzipleio[marzipleio$combinedoutcome %in% marziresults_plot$combinedoutcome,]
tableS4_p <- PACEpleio[PACEpleio$combinedoutcome %in% PACEresults_plot$combinedoutcome,]

# remove missings (as these do not have multiple snps anyways)
tableS4_7 <- tableS4_7[!is.na(tableS4_7$egger_intercept),]
tableS4_15 <- tableS4_15[!is.na(tableS4_15$egger_intercept),]
tableS4_m <- tableS4_m[!is.na(tableS4_m$egger_intercept),]
tableS4_p <- tableS4_p[!is.na(tableS4_p$egger_intercept),]

# bind
tableS4_7$age <- "Childhood DNAm"
tableS4_15$age <- "Adolescence DNAm"
tableS4_m$age <- "Young adulthood DNAm (triangulation with Marzi et al. (2019))"
tableS4_p$age  <- "Birth DNAm (triangulation with Ruehlmann et al. (2023))"
tableS4 <- rbind(tableS4_7, tableS4_15)
tableS4 <- rbind(tableS4, tableS4_m)
tableS4 <- rbind(tableS4, tableS4_p)

# select and clean up
tableS4 <- tableS4[,c("age", "exposure", "outcome", "egger_intercept", "se","pval")]
tableS4$egger_intercept <-  as.character(paste0(sprintf("%.2f", round(tableS4$egger_intercept,2))))
tableS4$se <-  as.character(paste0(sprintf("%.2f", round(tableS4$se,2))))
tableS4$pval <-  as.character(paste0(sprintf("%.3f", round(tableS4$pval,3))))

# write out
write.xlsx(tableS4, "YOURPATH/Analysis/FiguresTables/TableS7-2023-04-24.xlsx")

## ----------------------------------------------------------------------------------------
### TABLE S8

# retrieve significant CpGs
age7het$combinedoutcome <- paste0(age7het$outcome, " (", age7het$id.exposure, ")")
tableS5_7 <- age7het[age7het$combinedoutcome %in% age7results_plot$combinedoutcome,]
age15het$combinedoutcome <- paste0(age15het$outcome, " (", age15het$id.exposure, ")")
tableS5_15 <- age15het[age15het$combinedoutcome %in% age15results_plot$combinedoutcome,]
marzihet$combinedoutcome <- paste0(marzihet$outcome, " (", marzihet$id.exposure, ")")
tableS5_m <- marzihet[marzihet$combinedoutcome %in% marziresults_plot$combinedoutcome,]
PACEhet$combinedoutcome <- paste0(PACEhet$outcome, " (", PACEhet$id.exposure, ")")
tableS5_p <- PACEhet[PACEhet$combinedoutcome %in% PACEresults_plot$combinedoutcome,]

# bind
tableS5_7$age <- "Childhood DNAm"
tableS5_15$age <- "Adolescence DNAm"
tableS5_m$age <- "Young adulthood DNAm (triangulation with Marzi et al. (2019))"
tableS5_p$age  <- "Birth DNAm (triangulation with Ruehlmann et al. (2023))"
tableS5 <- rbind(tableS5_7, tableS5_15)
tableS5 <- rbind(tableS5, tableS5_m)
tableS5 <- rbind(tableS5, tableS5_p)

# select and clean up
tableS5 <- tableS5[tableS5$method != "MR Egger",]
tableS5 <- tableS5[,c("age", "exposure", "outcome","method", "Q", "Q_df","Q_pval")]
tableS5$Q <-  as.character(paste0(sprintf("%.2f", round(tableS5$Q,2))))
tableS5$Q_pval <-  as.character(paste0(sprintf("%.3f", round(tableS5$Q_pval,3))))

# write out
write.xlsx(tableS5, "YOURPATH/Analysis/FiguresTables/TableS8-2023-04-24.xlsx")

## ----------------------------------------------------------------------------------------
## Figure S1

# filter the ones with >SNP that were significant
age7Leave1out$combinedoutcome <- paste0(age7Leave1out$id.outcome, " (", age7Leave1out$id.exposure, ")")
figs1_7 <- age7Leave1out[age7Leave1out$combinedoutcome %in% age7results_plot$combinedoutcome[age7results_plot$nsnp > 2],]
age15Leave1out$combinedoutcome <- paste0(age15Leave1out$id.outcome, " (", age15Leave1out$id.exposure, ")")
figs1_15 <- age15Leave1out[age15Leave1out$combinedoutcome %in% age15results_plot$combinedoutcome[age15results_plot$nsnp > 2],]
marziLeave1out$combinedoutcome <- paste0(marziLeave1out$id.outcome, " (", marziLeave1out$id.exposure, ")")
figs1_m <- marziLeave1out[marziLeave1out$combinedoutcome %in% marziresults_plot$combinedoutcome[marziresults_plot$nsnp > 2],]
PACEleave1out$combinedoutcome <- paste0(PACEleave1out$id.outcome, " (", PACEleave1out$id.exposure, ")")
figs1_p <- PACEleave1out[PACEleave1out$combinedoutcome %in% PACEresults_plot$combinedoutcome[PACEresults_plot$nsnp > 2],]

# combine
leave1out <- rbind(figs1_7,figs1_15)
leave1out <- rbind(leave1out, figs1_m)
leave1out <- rbind(leave1out, figs1_p)

# fix names
leave1out[leave1out$id.outcome == "CKD","id.outcome"] <- "Chronic kidney disease"
leave1out[leave1out$id.outcome == "COPD","id.outcome"] <- "Chronic pulmonary disease"
leave1out[leave1out$id.outcome == "stroke","id.outcome"] <- "Stroke"
leave1out[leave1out$id.outcome == "scz","id.outcome"] <- "Schizophrenia"
leave1out[leave1out$id.outcome == "tourette","id.outcome"] <- "Tourette"
leave1out[leave1out$id.outcome == "ocd","id.outcome"] <- "Obsessive-compulsive disorder"
leave1out[leave1out$id.outcome == "anxiety","id.outcome"] <- "Anxiety disorders"
leave1out[leave1out$id.outcome == "exercise","id.outcome"] <- "Physical inactivity"
leave1out[leave1out$id.outcome == "CAD2022","id.outcome"] <- "Coronary artery disease"
leave1out$outcome <- leave1out$id.outcome

# plots
library(TwoSampleMR)
plots <- mr_leaveoneout_plot(leave1out)
plots[[1]] # to test

# write out
tiff("YOURPATH/Analysis/FiguresTables/FigureS3-2023-05-04.tiff", width = 6000, height = 3000, res = 400, pointsize = 4, type = "cairo")
plot_grid(plots[[1]],plots[[2]], plots[[3]], plots[[4]],
          plots[[5]],plots[[6]], plots[[7]], plots[[8]],
          plots[[9]],plots[[10]], plots[[11]],
          nrow = 3, ncol = 4)
dev.off()

