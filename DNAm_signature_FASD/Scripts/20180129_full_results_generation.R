

final2
load("~/NDN_FASD_I/svs2.RData")
load("~/NDN_FASD_I/svs.cauc.RData")
cluster1.samp <- read.table("~/NDN_FASD_I/Cluster1_Samples.txt")

mod <- model.matrix(~ as.factor(Group), data= pData(final2))
mod2 <- cbind(mod, svs2)
colnames(mod2) <- c("(Intercept)","FASD", paste("SV",1:15, sep =""))
head(mod2)

full.fit <- lmFit(final2, design = mod2)
full.eb <- eBayes(full.fit)
full.table <- topTable(full.eb, n=404030, adjust.method ="fdr", coef=2)
head(full.table)
full.dm <- full.table[(which(full.table$adj.P.Val <0.05)),] #1661


mod.clust <- model.matrix(~ as.factor(Group), data= pData(final2[,which(colnames(final2) %in% cluster1.samp$V1)]))
mod.clust2 <- cbind(mod.clust, svs.cauc)
dim(mod.clust2)

head(final2[,which(colnames(final2) %in% cluster1.samp$V1)])

cluster.fit <- lmFit(final2[,which(colnames(final2) %in% cluster1.samp$V1)], design = mod.clust2)
cluster.eb <- eBayes(cluster.fit)
cluster.table <- topTable(cluster.eb, coef = 2, n=Inf, adjust.method ="fdr")
clust.dm <- cluster.table[(which(cluster.table$P.Value <0.01)),] #5242


dm <- rownames(full.dm[which(rownames(full.dm) %in% rownames(clust.dm)),])
length(dm) #658



FASDsamples <- sampleNames(final2[,pData(final2)$Group=="FASD"])
CONTsamples <- sampleNames(final2[,pData(final2)$Group=="Control"])


library(lumi)
M.adj <- exprs(final2) - full.fit$coef[,3:17]%*%t(mod2[,3:17])
colnames(M.adj) <- sampleNames(final2)
B.adj <- m2beta(M.adj)
mean <- apply(B.adj, 1, mean, na.rm=T)
meanFASD <- apply(B.adj[,FASDsamples], 1, mean, na.rm=T)
meanCont <- apply(B.adj[,CONTsamples], 1, mean, na.rm=T)
deltas <- round(data.frame(mean.correctedBeta.FASD = meanFASD, mean.correctedBeta.Cont = meanCont,
                     delta.correctedBeta = meanFASD-meanCont),3)

deltas[which(rownames(deltas) == "cg10482532"),]
deltas.order <- deltas[order(rownames(deltas)),]
dim(full.order)

full.order <- full.table[order(rownames(full.table)),]
head(full.order)
identical(rownames(full.order), rownames(deltas.order))

clust.order <- cluster.table[order(rownames(cluster.table)),]
identical(rownames(full.order), rownames(clust.order))

results.full.table <- data.frame(Probe = rownames(full.order), P.Value = full.order$P.Value, fdr = full.order$adj.P.Val, 
                                 deltas.order,
                                 P.Value.Cluster1 = clust.order$P.Value)
head(results.full.table)
write.table(results.full.table, "~/NDN_FASD_I/full.results.table.txt", quote=F, row.names = F)

