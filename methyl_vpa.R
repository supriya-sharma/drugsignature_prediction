library(sva)

firstRows <- read.table("RReadyFinal.txt", header = TRUE, nrows = 5)
classes <- sapply(firstRows, class)
sFile <- read.table("RReadyFinal.txt", header = TRUE, colClasses = classes, nrow = 385000)
met <- sFile[complete.cases(sFile),]

meta <- read.csv("I157_Sample_Beadchip_Layout.csv")
sampleName <- paste(meta$Beadchip, meta$Strip, sep="_")
sampleID <- gsub(" ", "_", meta$Sample.ID)
searchTable <- cbind(sampleName, sampleID)
cellline <- as.character(sapply(names(met)[5:76],function(x){substr(x,2,nchar(x))}))

searchTable2 <- searchTable[match(cellline,searchTable[,1]),]

names(met)[5:76] <- searchTable2[,2]
idx <- paste(met[,1],met[,2],met[,3],met[,4],sep="_")

met2 <- met[,-(1:4)]
rownames(met2) <- idx

#reorder the samples in met2
met3 <- met2[,c(1:24,61:72,25:60,73:319)]

############
#PCA

pca <- prcomp(t(met[,-(1:4)]), retx=T, center=T, scale=T)
pc <- pca$x
col <- c(rep(1,24), rep(2,36), rep(1,12), rep(4,247))
col <- c(rep(1:2,length=12), rep(1,12), rep(3,36), rep(2,12),rep(4,247))

pdf("pca_methylation_noComBat.pdf")
plot(pc[,1], pc[,2], col=col,main="PCA for methylation data",pch=19)
legend("topright",legend=c("Mouse 110p5","Mouse 320p2","Cell lines","TCGA"),col=unique(col),pch=19)
dev.off()


pca <- prcomp(t(met3), retx=T, center=T, scale=T)
pc <- pca$x
col <- c(rep(1,26), rep(2,36), rep(3,247))
col <- c(rep(1:2,length=12), rep(1,12), rep(2,12), rep(3,36),rep(4,247))

pdf("pca_methylation_noComBat.pdf")
plot(pc[,1], pc[,2], col=col,main="PCA for methylation data",pch=19)
legend("topright",legend=c("Mouse 110p5","Mouse 320p2","Cell lines","TCGA"),col=unique(col),pch=19)
dev.off()

###########
# combat
batch <- c(rep("cellLine",36),rep("mouse",36), rep("tcga",247))
met_combat <- ComBat(dat=met3, batch, mod=NULL,numCovs=NULL, par.prior=TRUE,prior.plots=FALSE)
write.table(met3,file="RReadyFinal_v2.txt",quote=F,sep="\t")
write.table(met_combat,file="RReadyFinal_v2_combat.txt",quote=F,sep="\t")

## identify VPA 2h signature and VPA 6 h signature.
pb_cells <- met3[,37:72] 
vpa <- names(pb_cells)[grep("VPA", names(pb_cells))]
ctr <- names(pb_cells)[grep("ont", names(pb_cells))]

vpa_cell_ID <- c(ctr[c(1,3,5,7,9,11)],vpa[c(7:12,1:6)])
vpa_cells <- as.matrix(pb_cells[,vpa_cell_ID])

rownames(vpa_cells) <- idx

vpa_cells_logit <- log2((vpa_cells+0.001)/(1-(vpa_cells+0.001)))

library(limma)
treatment <- as.factor(rep(c("ctr","vpa2","vpa6"),each=6))
strain <- as.factor(rep(1:6,times=3))

design <-  model.matrix(~treatment+strain)
# fit <- lmFit(vpa_cells,design)
# fit <- eBayes(fit)
# topGenes_vpa2h <-topTable(fit,coef=2,number=200)
# topGenes_vpa6h <- topTable(fit, coef=3,number=200)

fit2 <- lmFit(vpa_cells_logit,design)
fit2 <- eBayes(fit2)

nTop <- 200
topGenes_vpa2h_2 <-topTable(fit2,coef=2,number=nTop)
topGenes_vpa6h_2 <- topTable(fit2, coef=3,number=nTop)

# associate methylation sites with genes
vpa2h <- topGenes_vpa2h_2[topGenes_vpa2h_2[,5]<0.05,]
vpa_id <- rownames(vpa2h)
vpa_gene <- NULL

for (i in 1:length(vpa_id)){
  vpa_gene <- c(vpa_gene, strsplit(vpa_id[i],split="_")[[1]][4])
}
vpa_gene_uniq <- unique(vpa_gene)
###########
### ASSIGN

library(ASSIGN, "/usr2/faculty/wej/R/x86_64-unknown-linux-gnu-library/2.15")
##VPA_2h
geneList_vpa2h <- rownames(topGenes_vpa2h_2)
S_matrix <- -fit2$coefficients[geneList_vpa2h,2]
B_vector <- fit2$coefficients[geneList_vpa2h,1]+fit2$coefficients[geneList_vpa2h,2]
Pi_matrix <- rep(0.95,nTop)

#TCGA
testData_sub <- met3[geneList_vpa2h,73:319]
testData_sub_logit <- log2((testData_sub+0.001)/(1-(testData_sub+0.001)))

#xenograph
testData_sub <- met3[geneList_vpa2h,1:36]
testData_sub_logit <- log2((testData_sub+0.001)/(1-(testData_sub+0.001)))

#110p5
xeno_110p5 <- testData_sub[,c(seq(1,12,by=2),13:24)]
#xeno_110p5_vpa <- xeno_110p5[,c(grep("ck",names(xeno_110p5)),grep("vpa",names(xeno_110p5)))][,c(1,2,3,7,8,9)]
testData_sub_logit <- log2((xeno_110p5+0.001)/(1-(xeno_110p5+0.001)))
testData_sub_logit <- testData_sub_logit[,order(names(testData_sub_logit))]                          

#320p2
xeno_320p2 <- testData_sub[,c(seq(2,12,by=2),25:36)]
#xeno_320p2_vpa <- xeno_320p2[,c(grep("ck",names(xeno_320p2)),grep("vpa",names(xeno_320p2)))][,c(1,2,3,7,8,9)]

testData_sub_logit <- log2((xeno_320p2+0.001)/(1-(xeno_320p2+0.001)))
testData_sub_logit <- testData_sub_logit[,order(names(testData_sub_logit))]



##VPA_6h
geneList_vpa6h <- rownames(topGenes_vpa6h_2)
testData_sub <- met[geneList_vpa6h,77:323]
testData_sub_logit <- log2((testData_sub+0.001)/(1-(testData_sub+0.001)))

S_matrix <- -fit2$coefficients[geneList_vpa6h,3]
B_vector <- fit2$coefficients[geneList_vpa6h,1]+fit2$coefficients[geneList_vpa6h,3]
Pi_matrix <- rep(0.95,nTop)

#TCGA
testData_sub <- met3[geneList_vpa6h,73:319]
testData_sub_logit <- log2((testData_sub+0.001)/(1-(testData_sub+0.001)))

#xenograph
testData_sub <- met3[geneList_vpa6h,1:36]
testData_sub_logit <- log2((testData_sub+0.001)/(1-(testData_sub+0.001)))

#110p5
xeno_110p5 <- testData_sub[,c(seq(1,12,by=2),13:24)]
#xeno_110p5_vpa <- xeno_110p5[,c(grep("ck",names(xeno_110p5)),grep("vpa",names(xeno_110p5)))][,c(1,2,3,7,8,9)]
testData_sub_logit <- log2((xeno_110p5+0.001)/(1-(xeno_110p5+0.001)))
testData_sub_logit <- testData_sub_logit[,order(names(testData_sub_logit))]                          

#320p2
xeno_320p2 <- testData_sub[,c(seq(2,12,by=2),25:36)]
#xeno_320p2_vpa <- xeno_320p2[,c(grep("ck",names(xeno_320p2)),grep("vpa",names(xeno_320p2)))][,c(1,2,3,7,8,9)]

testData_sub_logit <- log2((xeno_320p2+0.001)/(1-(xeno_320p2+0.001)))
testData_sub_logit <- testData_sub_logit[,order(names(testData_sub_logit))]

#test1: adaptive_B=F, adaptive_S=F, mixture_beta=F
mcmc.chain <- assign.mcmc(Y = testData_sub_logit, Bg = B_vector, X = S_matrix, Delta_prior_p = Pi_matrix, iter=2000, adaptive_B=F, adaptive_S=F, mixture_beta=F, p_beta = 0.5)
mcmc.pos.mean1 <- assign.summary(test=mcmc.chain, burn_in=1000, iter=2000, adaptive_B=F, adaptive_S=F,mixture_beta=F)

#test2: adaptive_B=T, adaptive_S=F, mixture_beta=F
mcmc.chain <- assign.mcmc(Y = testData_sub_logit, Bg = B_vector, X = S_matrix, Delta_prior_p = Pi_matrix, iter=2000, adaptive_B=T, adaptive_S=F, mixture_beta=F, p_beta = 0.5)
mcmc.pos.mean2 <- assign.summary(test=mcmc.chain, burn_in=1000, iter=2000, adaptive_B=T, adaptive_S=F,mixture_beta=F)

#test3: adaptive_B=T, adaptive_S=T, mixture_beta=F
mcmc.chain <- assign.mcmc(Y = testData_sub_logit, Bg = B_vector, X = S_matrix, Delta_prior_p = Pi_matrix, iter=2000, adaptive_B=T, adaptive_S=T, mixture_beta=F, p_beta = 0.5)
mcmc.pos.mean3 <- assign.summary(test=mcmc.chain, burn_in=1000, iter=2000, adaptive_B=T, adaptive_S=T,mixture_beta=F)

#test4: adaptive_B=T, adaptive_S=T, mixture_beta=T
mcmc.chain <- assign.mcmc(Y = testData_sub_logit, Bg = B_vector, X = S_matrix, Delta_prior_p = Pi_matrix, iter=2000, adaptive_B=T, adaptive_S=T, mixture_beta=T, p_beta = 0.5)
mcmc.pos.mean4 <- assign.summary(test=mcmc.chain, burn_in=1000, iter=2000, adaptive_B=T, adaptive_S=T,mixture_beta=T)
vpa_pa <- mcmc.pos.mean4$kappa_pos
row.names(vpa_pa) <- names(testData_sub)
write.csv(vpa_pa, file="methylation_tcga_vpa_all.csv")

# plots and tables
label <- as.factor(c(rep("normal",32),rep("tumor",215)))
label <- as.factor(c(rep("ck",6),rep("vpa_treated",6),rep("zeb_treated",6)))

pdf("methylation_tcga_vpa.pdf")
par(mfrow=c(2,2))
boxplot(mcmc.pos.mean1$beta ~ label,ylab="vpa signature",main="test1")
boxplot(mcmc.pos.mean2$beta ~ label,ylab="vpa signature",main="test2")
boxplot(mcmc.pos.mean3$beta ~ label,ylab="vpa signature",main="test3")
boxplot(mcmc.pos.mean4$kappa ~ label,ylab="vpa signature",main="test4")
dev.off()

pa <- matrix(mcmc.pos.mean4$kappa,nrow=ncol(testData_sub_logit),1);rownames(pa) <- names(testData_sub_logit)
pa1 <- pa[order(rownames(pa)),1,drop=F]
pa2 <- pa1[c(grep("11A",rownames(pa1))-1, grep("11B",rownames(pa1))-1,grep("11A",rownames(pa1)), grep("11B",rownames(pa1))),1,drop=F]
dim(pa2) <- c(32,2)
rownames(pa2) <- sapply(rownames(pa1)[1:32],function(x){substr(x,1,12)})
colnames(pa2) <- c("tumor", "normal")

write.csv(pa2,file="methylation_tcga_vpa.csv")

###
coeff <- mcmc.pos.mean4$kappa
rownames(coeff) <- names(testData_sub_logit)
colnames(coeff) <- "vpa_6h_signature"
write.csv(coeff,file="110p5_vpa6h_methylation_signature.csv")

coeff <- mcmc.pos.mean4$kappa
rownames(coeff) <- names(testData_sub_logit)
colnames(coeff) <- "vpa_6h_signature"
write.csv(coeff,file="320p2_vpa6h_methylation_signature.csv")

pdf("methylation_110p5_vpa_6h.pdf")
boxplot(mcmc.pos.mean4$kappa ~ label,ylab="vpa signature",main="110p5_vpa_6h_methylation_signature.pdf")
dev.off()                           

pdf("methylation_320p2_vpa_6h.pdf")
boxplot(mcmc.pos.mean4$kappa ~ label,ylab="vpa signature",main="320p2_vpa_6h_methylation_signature.pdf")
dev.off()                           

pdf("methylation_tcga_vpa_6h.pdf")
boxplot(mcmc.pos.mean4$kappa ~ label,ylab="vpa signature",main="TCGA_vpa_6h_methylation_signature.pdf")
dev.off()                           
                           