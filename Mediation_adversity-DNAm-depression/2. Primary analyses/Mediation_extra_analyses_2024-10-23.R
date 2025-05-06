


#*Brain-blood correlations of mediation CpGs*
library(ggplot2)
alspac <- readxl::read_xls("~/Dropbox (Partners HealthCare)/Dunn_shared/MGH-GenR-mediation/Dummy_data_2021-11-04/data/Examples/Mediation_sites_all_2021-06-04.xls")
head(alspac)
tail(alspac)
alspac <- alspac[-71,]

brainBlood <- data.table::fread("~/Dropbox (Partners HealthCare)/Dunn_shared/OneMind/SLCMA-Bristol-MGH/Analyses/450K_annotations/ImageCpG/Illumina_450K_covar")

brainBloodMed <- brainBlood[match(alspac$CpG,brainBlood$cgid), ]
mean(brainBloodMed$rho.br.b, na.rm=T)

sum(brainBloodMed$rho.br.b >0.8, na.rm=T)

brainBloodMedCor <- rbind(data.frame(alspac, 
                                         BB.rho = brainBloodMed$rho.br.b,
                                         BB.pval = brainBloodMed$p.rho.br.b,
                                         BB.analysis = "Brain-Blood"),
                              data.frame(alspac, 
                                         BB.rho = brainBloodMed$rho.br.s,
                                         BB.pval = brainBloodMed$p.rho.br.s,
                                         BB.analysis = "Brain-Saliva"),
                              data.frame(alspac, 
                                         BB.rho = brainBloodMed$rho.bl.s,
                                         BB.pval = brainBloodMed$p.rho.bl.s,
                                         BB.analysis = "Blood-Saliva"))


adv.labels <- c("Family instability",
                "Financial stress",
                "Mom psychopathology",
                "One adult household",
                "Parental cruelty")

brainBloodMedCor$Adversity <- factor(brainBloodMedCor$Adversity, 
                                     levels = c("parc","abuse","mompsych",
                                                "oneadult","r_faminst",
                                                "Fscore","nbhqual"))
brainBloodMedCor$BB.analysis <- factor(brainBloodMedCor$BB.analysis,
                                       levels = c("Brain-Blood","Brain-Saliva","Blood-Saliva"))
table(brainBloodMedCor$Adversity)/3
ggplot(brainBloodMedCor, aes(x = Adversity, y= BB.rho, fill=BB.analysis))+
  #geom_jitter(aes(col = BB.analysis, group = BB.analysis))
  #geom_violin(na.rm = T, width= 1, alpha =0.7)+
  geom_boxplot(na.rm = T, width= 0.5, alpha =0.7, position='dodge')+
  #geom_dotplot(binaxis = "y", stackdir = "center", dotsize = 0.5, position = "dodge")+
  theme_classic()+
  scale_fill_viridis(discrete=T, "Tissues")+
  geom_hline(yintercept = 0, linetype =2)+
  xlab("Adversity")+
  ylab("Correlation coefficient (rho)")+
  theme(panel.grid.major.y = element_line(linetype=3),
        axis.text.x = element_text(size =12, angle=45, hjust=1, color="black"),
        axis.text.y = element_text(size =10, color="black"),
        axis.title = element_text(size =14, color="black"),
        legend.text = element_text(size =10, color="black"),
        legend.title = element_text(size =12, color="black"),
        strip.text = element_text(size=12, color="black"))+
  scale_x_discrete(labels = c("Caregiver physical\nor emotional abuse\n(5 CpGs)", 
                               "Sexual or physical\nabuse (by anyone)\n(6 CpGs)", 
                               "Maternal\npsychopathology\n(9 CpGs)", 
                               "One adult in\nthe household\n(10 CpGs)", 
                               "Family\ninstability\n(8 CpGs)", 
                               "Financial\nhardship\n(17 CpGs)", 
                               "Neighborhood\ndisadvantage\n(15 CpGs)"))


ggplot(brainBloodMedCor, aes(y=BB.rho, x = BB.analysis, fill = BB.analysis))+
  geom_hline(yintercept = 0, linetype =3)+
  geom_boxplot(na.rm=T, alpha = 0.7)+
  theme_classic()+
  scale_fill_viridis(discrete=T, "Tissues")+
  ylab("Correlation coefficient (rho)")+
  xlab("Compared tissues")+
  theme(panel.grid.major.y = element_line(linetype=3),
        axis.text.x = element_text(size =12, angle=45, hjust=1),
        axis.text.y = element_text(size =10),
        axis.title = element_text(size =14),
        legend.text = element_text(size =10),
        legend.title = element_text(size =12),
        strip.text = element_text(size=12))

grep("cg27423959", alspac$CpG)

## trying new figure
load("~/Dropbox (Partners HealthCare)/Dunn_shared/ALSPAC/Epigenetic/Causal Inference/Mediation Analysis/output/Mediation_sites_2022-11-28.rdata")
View(med.data)
library(tidyr)
library(ggplot2)
library(dplyr)
alspac <- med.data
alspac$type <- ifelse(alspac$alpha>0 & alspac$beta>0,"pos.pos",
                      ifelse(alspac$alpha<0 & alspac$beta<0,"neg.neg",
                             ifelse(alspac$alpha>0 & alspac$beta<0,"pos.neg",
                                    ifelse(alspac$alpha<0 & alspac$beta>0,"neg.pos",NA))))
alspac <- alspac[-71,]

table(alspac$alpha<0, alspac$beta<0)


colnames(alspac)
alspac$type <- factor(alspac$type, levels =c("pos.pos","neg.neg","pos.neg","neg.pos"))
alspac$type2 <- factor(alspac$type, levels = c("Risk-increasing", "pos.pos","neg.neg",
                                               "Risk-suppressing","pos.neg","neg.pos"))
alspac$Adversity <- factor(alspac$Adversity,
                           levels = c("parc","abuse","mompsych","nbhqual","oneadult","r_faminst","Fscore"))


advNames <- c("Caregiver physical\nor emotional abuse","Sexual/physical\nabuse (by anyone)","Maternal\npsychopathology",
  "Neighborhood\ndisadvantage","One-adult\nhousehold","Family\ninstability",
  "Financial\nhardship")
names(advNames) <- levels(alspac$Adversity)
alspac$Adversity2 <- advNames[names(advNames)[alspac$Adversity]]
alspac$Adversity2 <- factor(alspac$Adversity2, levels = advNames)

alspac$SP <- ifelse(alspac$Adversity %in% c("parc","abuse","mompsych","nbhqual"), 
                    "Very early childhood","Early childhood")
alspac$SP <- factor(alspac$SP, levels =c("Very early childhood","Early childhood"))
alspac <- alspac[order(abs(alspac$alphaxbeta/alspac$`sum_delta + alpha*beta`), decreasing = T),]
alspac$Adversity3 <- factor(alspac$Adversity2, levels = rev(levels(alspac$Adversity2)))
alspac$SP <- factor(alspac$SP,levels = c("Early childhood","Very early childhood"))
netEffect <- data.frame(Adversity3 = levels(alspac$Adversity2),
                        medEffect = c("-27%","-73%","-22%","-13%","-71%","+10%","+30%"),
                        SP = c(rep("Very early childhood",4), rep("Early childhood",3)))


# ggplot(alspac, aes(y= Adversity3, x= (alphaxbeta/`sum_delta + alpha*beta`)*100, fill =type2, col=type2))+
#   geom_vline(xintercept = c(25, 50, 75, 100, -25, -50,-75, -100), col="grey", size =0.1)+
#   geom_bar(stat='identity', alpha=0.77, width=0.77)+
#   geom_vline(xintercept = 0)+
#   theme_classic()+
#   facet_grid(rows = "SP",scales = "free_y", space="free")+
#   scale_fill_manual("Mediator pattern", values = c(NA,"#ea4040", "#b70d0d",NA,"#66c1e9","#0088c5"),
#                     labels = c(expression(bold("Risk-increasing")),"Positive-Positive (n=19)","Negative-Negative (n=12)",
#                                expression(bold("Risk-suppressing")),"Positive-Negative (n=23)","Negative-Positive (n=16)"),
#                     na.value="white", drop=F)+
#   scale_color_manual("Mediator pattern", values = c(NA,"#393939", "#393939",NA,"#393939","#393939"),
#                     labels = c(expression(bold("Risk-increasing")),"Positive-Positive (n=19)","Negative-Negative (n=12)",
#                                expression(bold("Risk-suppressing")),"Positive-Negative (n=23)","Negative-Positive (n=16)"),
#                     na.value="white", drop=F)+
#   theme(axis.text = element_text(size = 10),
#         axis.title = element_text(size = 12),
#         legend.text = element_text(size = 10),
#         legend.title = element_text(size = 12),
#         legend.text.align = 0.5,
#         legend.title.align = 0.7)+
#   xlab("Percent mediated effect")+
#   ylab("")+
#   xlim(-115, 115)+
#   #geom_text(aes(label = round(sum((alphaxbeta/`sum_delta + alpha*beta`)*100))))
#   geom_text(data = netEffect, aes(y = Adversity3, label= medEffect, x = 115), inherit.aes = F)

#flipped - figure for paper
alspac$SP <- factor(alspac$SP,levels = c("Very early childhood", "Early childhood"))
netEffect$SP <- factor(netEffect$SP, levels = c("Very early childhood", "Early childhood"))
alspac$Adversity3

p <- ggplot(alspac, aes(x= Adversity2, y= (alphaxbeta/`sum_delta + alpha*beta`)*100, fill =type2, col=type2))+
  geom_hline(yintercept = c(25, 50, 75, 100, -25, -50,-75, -100), col="grey", size =0.1)+
  geom_bar(stat='identity', alpha=0.77, width=0.77)+
  geom_hline(yintercept = 0)+
  theme_classic()+
  facet_grid(~SP,scales = "free_x", space="free")+
  #facet_wrap(~SP, scales = "free_x")+
  scale_fill_manual("Mediator pattern", values = c(NA,"#ea4040", "#b70d0d",NA,"#66c1e9","#0088c5"),
                    labels = c(expression(bold("Risk-increasing")),"Positive-Positive (n=19)","Negative-Negative (n=12)",
                               expression(bold("Protective")),"Positive-Negative (n=23)","Negative-Positive (n=16)"),
                    na.value="white", drop=F)+
  scale_color_manual("Mediator pattern", values = c(NA,"#393939", "#393939",NA,"#393939","#393939"),
                     labels = c(expression(bold("Risk-increasing")),"Positive-Positive (n=19)","Negative-Negative (n=12)",
                                expression(bold("Protective")),"Positive-Negative (n=23)","Negative-Positive (n=16)"),
                     na.value="white", drop=F)+
  theme(axis.text = element_text(size = 10),
        axis.text.x = element_text(angle=45, hjust=1),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.text.align = 0.5,
        legend.title.align = 0.7)+
  ylab("Indirect (mediated) effect")+
  xlab("")+
  ylim(-115, 115)+
  #geom_text(aes(label = round(sum((alphaxbeta/`sum_delta + alpha*beta`)*100))))
  geom_text(data = netEffect, aes(x = Adversity3, label= medEffect, y = 115), inherit.aes = F)
p

ggsave(p, file = "~/Dropbox (Partners HealthCare)/Dunn_shared/ALSPAC/Epigenetic/Paper Draft - mediation/Proofs/Figure1_2024-10-23.eps", 
       dpi = 300, width=10, height=6, device=cairo_ps)

write.csv(alspac[,!colnames(alspac) %in% c("Adversity2", "Adversity3")], 
          file="~/Dropbox (Partners HealthCare)/Dunn_shared/ALSPAC/Epigenetic/Paper Draft - mediation/Submission 5 - Nature Mental Health editorial/Data/Mediation_figure1_data_2024-08-27.csv",
          quote=F, row.names = F)

#no panels
# alspac$Adversity4 <- factor(alspac$Adversity3, 
#                             levels =netEffect$Adversity3[order(as.numeric(gsub("%","",netEffect$medEffect)))])
# 
# ggplot(alspac, aes(y= Adversity3, x= (alphaxbeta/`sum_delta + alpha*beta`)*100, fill =type2, col=type2))+
#   geom_vline(xintercept = c(25, 50, 75, 100, -25, -50,-75, -100), col="grey", size =0.1)+
#   geom_bar(stat='identity', alpha=0.77, width=0.77)+
#   geom_vline(xintercept = 0)+
#   theme_classic()+
#   #facet_grid(rows = "SP",scales = "free_y", space="free")+
#   scale_fill_manual("Mediator pattern", values = c(NA,"#ea4040", "#b70d0d",NA,"#66c1e9","#0088c5"),
#                     labels = c(expression(bold("Risk-increasing")),"Positive-Positive (n=19)","Negative-Negative (n=12)",
#                                expression(bold("Risk-suppressing")),"Positive-Negative (n=23)","Negative-Positive (n=16)"),
#                     na.value="white", drop=F)+
#   scale_color_manual("Mediator pattern", values = c(NA,"#393939", "#393939",NA,"#393939","#393939"),
#                      labels = c(expression(bold("Risk-increasing")),"Positive-Positive (n=19)","Negative-Negative (n=12)",
#                                 expression(bold("Risk-suppressing")),"Positive-Negative (n=23)","Negative-Positive (n=16)"),
#                      na.value="white", drop=F)+
#   theme(axis.text = element_text(size = 10),
#         axis.title = element_text(size = 12),
#         legend.text = element_text(size = 10),
#         legend.title = element_text(size = 12),
#         legend.text.align = 0.5,
#         legend.title.align = 0.7)+
#   xlab("Percent mediated effect")+
#   ylab("")+
#   xlim(-115, 115)+
#   #geom_text(aes(label = round(sum((alphaxbeta/`sum_delta + alpha*beta`)*100))))
#   geom_text(data = netEffect, aes(y = Adversity3, label= medEffect, x = 115), inherit.aes = F)
# 
# 
# 
# 
# alspac %>% group_by(Adversity) %>% summarize(sum = round(sum((alphaxbeta/`sum_delta + alpha*beta`)*100)))
# alspac$`Sum alpha*beta`/alspac$delta

# library(wesanderson)
# ggplot(alspac, aes(x= Adversity, fill =type))+
#   geom_bar(stat='count', alpha=0.8, col='black')+
#   theme_classic()+
#   scale_fill_manual(values = c("darkred", "blue","darkblue","red"))+
#   geom_hline(yintercept = 0)
# 
# 
# t <- alspac[alspac$Adversity=="r_faminst",]
# sum(t$alphaxbeta)
# 
# alspac%>% group_by(Adversity) %>% summarize(sum = sum(alphaxbeta),
#                                             abs.sum = sum(Estimated.mediation.effect..abs.val.of.alphaxbeta.),
#                                             sum.alphabeta = mean(Sum.alpha.beta),
#                                             total = mean(sum_delta...alpha.beta))


#sample size FFCWS
load("/Users/aq868/Partners HealthCare Dropbox/Alexandre Lussier/Dunn_shared/ALSPAC/Epigenetic/Paper Draft - mediation/Submission 4 - Nature Mental Health R&R/Analyses/SLCMA_replication/FFCWS_site_replication(array)_2024-03-25.Rdata", verbose=T)
allMediation[allMediation$CpG =="cg09191574",]
t <- table(allMediation$n, allMediation$Adversity)
t[,-grep(".ever", colnames(t))]
table(allMediation$CpG[allMediation$n==773])

temp <- allMediation[allMediation$Mediation=="Total Effect",]
temp <- temp[!temp$Adversity =="oneadult5",]
temp <- temp[-grep(".ever",temp$Adversity),]
temp[temp$Adversity =="par3",]
unique(temp$Adversity)

temp[temp$Adversity =="par3" & temp$CpG %in% alspac$CpG[alspac$Adversity=="parc"],]
temp[temp$Adversity =="dep3" & temp$CpG %in% alspac$CpG[alspac$Adversity=="mompsych"],]
temp[temp$Adversity =="oneadult3" & temp$CpG %in% alspac$CpG[alspac$Adversity=="oneadult"],]
temp[temp$Adversity =="faminst5" & temp$CpG %in% alspac$CpG[alspac$Adversity=="r_faminst"],]
temp[temp$Adversity =="fh5" & temp$CpG %in% alspac$CpG[alspac$Adversity=="Fscore"],]


