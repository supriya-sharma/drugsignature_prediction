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
##col <- c(rep(1,24), rep(2,36), rep(1,12), rep(4,247))
##col <- c(rep(1:2,length=12), rep(1,12), rep(3,36), rep(2,12),rep(4,247))


col.110p5_ck <- which(substr(rownames(pc),1,8) == "110p5_ck")
col.110p5_zeb <- which(substr(rownames(pc),1,9) == "110p5_zeb")
col.110p5_vpa <- which(substr(rownames(pc),1,9) == "110p5_vpa")

col.320p2_ck <- which(substr(rownames(pc),1,8) == "320p2_ck")
col.320p2_zeb <- which(substr(rownames(pc),1,9) == "320p2_zeb")
col.320p2_vpa <- which(substr(rownames(pc),1,9) == "320p2_vpa")

col.TCGA <- which(substr(rownames(pc),1,4) == "TCGA")

col.cell.lines <- (1:nrow(pc))[-c(col.110p5_ck, col.110p5_zeb, col.110p5_vpa,
                                  col.320p2_ck, col.320p2_zeb, col.320p2_vpa,
                                  col.TCGA)]

#col <- c(rep(1:2,length=12), rep(1,12), rep(3,36), rep(2,12),rep(4,247))
my.col <- c(rep(1, length(col.110p5_ck)), 
            rep(2, length(col.110p5_zeb)),
            rep(3, length(col.110p5_vpa)),
            rep(4, length(col.320p2_ck)),
            rep(5, length(col.320p2_zeb)),
            rep(6, length(col.320p2_vpa)),
            rep(7, length(col.TCGA)),
            rep(8, length(col.cell.lines))
            
)
pdf("pca_methylation$$_noComBat.pdf")
plot(pc[c(col.110p5_ck, 
          col.110p5_zeb, 
          col.110p5_vpa,
          col.320p2_ck, 
          col.320p2_zeb, 
          col.320p2_vpa,
          col.TCGA, 
          col.cell.lines),1], 
     pc[c(col.110p5_ck, 
          col.110p5_zeb, 
          col.110p5_vpa,
          col.320p2_ck, 
          col.320p2_zeb, 
          col.320p2_vpa,
          col.TCGA, 
          col.cell.lines),2], 
     col=my.col,main="PCA for methylation data",pch=19, 
     xlab="", ylab="")

#pdf("pca_methylation$$_noComBat.pdf")
legend("topright",legend=c("110p5ck","110p5zeb","110p5vpa","320p2ck","320p2zeb","320p2vpa","TCGA","Cell lines),col=unique(col),pch=19)
dev.off()

)
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
## With batch adjusted data
pb_cells <- met_combat[,37:72] 
vpa <- names(pb_cells)[grep("VPA", names(pb_cells))]
ctr <- names(pb_cells)[grep("ont", names(pb_cells))]

vpa_cell_ID <- c(ctr[c(1,3,5,7,9,11)],vpa[c(7:12,1:6)])
vpa_cells <- as.matrix(pb_cells[,vpa_cell_ID])

rownames(vpa_cells) <- idx
vpa_cells[vpa_cells<0] = 0
vpa_cells[vpa_cells>.998] = .998

vpa_cells_logit <- log2((vpa_cells+0.001)/(1-(vpa_cells+0.001)))

library(limma)
treatment <- as.factor(rep(c("ctr","vpa2","vpa6"),each=6))
#treatment <- as.factor(rep(c("ctr","vpa","vpa"),each=6))
strain <- as.factor(rep(1:6,times=3))

design <-  model.matrix(~treatment+strain)
# fit <- lmFit(vpa_cells,design)
# fit <- eBayes(fit)
#topGenes_vpa2h <-topTable(fit,coef=2,number=200)
#topGenes_vpa6h <- topTable(fit, coef=3,number=200)

fit2 <- lmFit(vpa_cells_logit,design)
fit2 <- eBayes(fit2)

# 5000 gene selection
nTop <- 5000
topGenes_vpa2h_2 <-topTable(fit2,coef=2,number=nTop)
topGenes_vpa6h_2 <- topTable(fit2, coef=3,number=nTop)

## write.csv(topGenes_vpa6h_2, file="vpa_diff_Methyl_GeneList_5000.csv")
       
#  associate methylation sites with genes vpa2
vpa2h <- topGenes_vpa2h_2[topGenes_vpa2h_2[,5]<0.05,]
vpa2_id <- rownames(vpa2h)

vpa2_gene <- NULL
for (i in 1:length(vpa2_id)){
  vpa2_gene <- c(vpa2_gene, strsplit(vpa2_id[i],split="_")[[1]][4])
}
vpa2_gene_uniq <- unique(vpa2_gene)
keep2 = (vpa2_gene != "NONE") & (!duplicated(vpa2_gene))   ## Evan remove diuplicated probes and controle probes ("NONE")
topGenes_vpa2h_2_keep = topGenes_vpa2h_2[keep2,]

#write.csv(topGenes_vpa2h_2, file="VPA_diff_Me_Unique_GeneList_5000.csv")

# associate methylation sites with genes vpa6
vpa6h <- topGenes_vpa6h_2[topGenes_vpa6h_2[,5]<0.05,]
vpa6_id <- rownames(vpa6h)
vpa6_id <- vpa6h$ID
vpa6_gene <- NULL
for (i in 1:length(vpa6_id)){
  vpa6_gene <- c(vpa6_gene, strsplit(vpa6_id[i],split="_")[[1]][4])
}
vpa6_gene_uniq <- unique(vpa6_gene)
keep6 = (vpa6_gene != "NONE") & (!duplicated(vpa6_gene))   ## Evan remove diuplicated probes and controle probes ("NONE")
topGenes_vpa6h_2_keep = topGenes_vpa6h_2[keep6,][1:5000,]

write.csv(topGenes_vpa6h_2_keep, file="VPA_diff_Me_Unique_GeneList_5000.csv")

#### read in gene list
### reduce topGenes_

######################################

##### Evan's edits#################
### ASSIGN

library(ASSIGN, "/usr2/faculty/wej/R/x86_64-unknown-linux-gnu-library/2.15")


# 5000 gene selection
nTop <- 5000

topGenes_vpa2h_2 <-topTable(fit2,coef=2,number=nTop)
topGenes_vpa6h_2 <- topTable(fit2, coef=3,number=nTop)
                           
#  associate methylation sites with genes vpa2
vpa2h <- topGenes_vpa2h_2[topGenes_vpa2h_2[,5]<0.05,]

vpa2_id <- rownames(vpa2h)

vpa2_gene <- NULL
for (i in 1:length(vpa2_id)){
vpa2_gene <- c(vpa2_gene, strsplit(vpa2_id[i],split="_")[[1]][4])
}
vpa2_gene_uniq <- unique(vpa2_gene)
keep2 = (vpa2_gene != "NONE") & (!duplicated(vpa2_gene))   ## Evan remove diuplicated probes and controle probes ("NONE")
                           
topGenes_vpa2h_2_keep = topGenes_vpa2h_2[keep2,]  

#write.csv(topGenes_vpa2h_2, file="VPA_diff_Me_Unique_GeneList_5000.csv")
                           
# associate methylation sites with genes vpa6
vpa6h <- topGenes_vpa6h_2[topGenes_vpa6h_2[,5]<0.05,]
                          
vpa6_id <- rownames(vpa6h)

vpa6_gene <- NULL
for (i in 1:length(vpa6_id)){
vpa6_gene <- c(vpa6_gene, strsplit(vpa6_id[i],split="_")[[1]][4])
}
vpa6_gene_uniq <- unique(vpa6_gene)
keep6 = (vpa6_gene != "NONE") & (!duplicated(vpa6_gene))   ## Evan remove diuplicated probes and controle probes ("NONE")
                           
topGenes_vpa6h_2_keep = topGenes_vpa6h_2[keep6,] 

write.csv(topGenes_vpa6h_2_keep, file="VPA_diff_Me_Unique_GeneList_5000.csv")
                           

topGenes_unique <- read.csv("correlation_unique_500_geneList.csv", header = TRUE)
                           
### Link exp/meth neg corr genes with limma result above:
#get the probe names (from topGenes_vpa6h_2_keep) for the genes in the correlated list

selected=NULL
for (i in rownames(topGenes_vpa6h_2_keep)){
selected <- c(selected, strsplit(i,split= "_")[[1]][4]%in%topGenes_unique$genes)
}
topCor=rownames(topGenes_vpa6h_2_keep)[selected]
                           
S_matrix <- -fit2$coefficients[rownames(fit2$coefficients)%in%topCor,2]
B_vector <- fit2$coefficients[rownames(fit2$coefficients)%in%topCor,1]+fit2$coefficients[rownames(fit2$coefficients)%in%topCor,2]
Pi_matrix <- rep(0.95,nrow(topGenes_unique))
                           
##VPA_6h
                           
geneList_vpa6h <- rownames(topGenes_vpa6h_2_keep)
                           
S_matrix <- -fit2$coefficients[rownames(fit2$coefficients)%in%topCor,3]
B_vector <- fit2$coefficients[rownames(fit2$coefficients)%in%topCor,1]+fit2$coefficients[rownames(fit2$coefficients)%in%topCor,3]
Pi_matrix <- rep(0.95,length(S_matrix))
                           
#TCGA
                         
testData_sub_TCGA <-met_combat[topCor,73:319]
testData_sub_TCGA[testData_sub_TCGA<0] = 0
testData_sub_TCGA[testData_sub_TCGA>.998] = .998
testData_sub_TCGA_logit <- log2((testData_sub_TCGA+0.001)/(1-(testData_sub_TCGA+0.001)))


#test4: adaptive_B=T, adaptive_S=T, mixture_beta=T
mcmc.chain <- assign.mcmc(Y = testData_sub_TCGA_logit , Bg = B_vector, X = S_matrix, Delta_prior_p = Pi_matrix, iter=2000, adaptive_B=F, adaptive_S=F, mixture_beta=T, p_beta = 0.5)
mcmc.pos.mean4 <- assign.summary(test=mcmc.chain, burn_in=1000, iter=2000, adaptive_B=F, adaptive_S=F,mixture_beta=T)
vpa_pa <- mcmc.pos.mean4$beta_pos
##############################################################################
#xenograph
testData_sub_xenograph <- met_combat[geneList_vpa2h,1:36]
testData_sub_xenograph_logit <- log2((testData_sub_xenograph+0.001)/(1-(testData_sub_xenograph+0.001)))

#110p5
xeno_110p5 <- testData_sub_xenograph[,c(seq(1,12,by=2),13:24)]
xeno_110p5_vpa <- xeno_110p5[,c(grep("ck",names(xeno_110p5)),grep("vpa",names(xeno_110p5)))][,c(1,2,3,7,8,9)]
testData_sub_110p5_logit <- log2((xeno_110p5+0.001)/(1-(xeno_110p5+0.001)))
testData_sub_110p5_logit <- testData_sub_110p5_logit[,order(names(testData_sub_110p5_logit))]                          

#320p2
##xeno_320p2 <- testData_sub_xenograph[,c(seq(2,12,by=2),25:36)]
#xeno_320p2_vpa <- xeno_320p2[,c(grep("ck",names(xeno_320p2)),grep("vpa",names(xeno_320p2)))][,c(1,2,3,7,8,9)]

##testData_sub_320p2_logit <- log2((xeno_320p2+0.001)/(1-(xeno_320p2+0.001)))
##testData_sub_320p2_logit <- testData_sub_320p2_logit[,order(names(testData_sub_320p2_logit))]



##VPA_6h
geneList_vpa6h <- rownames(topGenes_vpa6h_2)
testData_sub_vpa6h <- met_combat[geneList_vpa6h,77:323]
testData_sub_vpa6h_logit <- log2((testData_sub_vpa6h+0.001)/(1-(testData_sub_vpa6h+0.001)))

S_matrix <- -fit2$coefficients[geneList_vpa6h,3]
B_vector <- fit2$coefficients[geneList_vpa6h,1]+fit2$coefficients[geneList_vpa6h,3]
Pi_matrix <- rep(0.95,nTop)

#TCGA
testData_sub_TCGA <- met_combat[geneList_vpa6h,73:319]
testData_sub_TCGA_logit <- log2((testData_sub_TCGA+0.001)/(1-(testData_sub_TCGA+0.001)))

#xenograph
testData_sub_xenograph <- met_combat[geneList_vpa6h,1:36]
testData_sub_xenograph_logit <- log2((testData_sub_xenograph+0.001)/(1-(testData_sub_xenograph+0.001)))

#110p5
xeno_110p5 <- testData_sub_xenograph[,c(seq(1,12,by=2),13:24)]
xeno_110p5_vpa <- xeno_110p5[,c(grep("ck",names(xeno_110p5)),grep("vpa",names(xeno_110p5)))][,c(1,2,3,7,8,9)]
testData_sub_110p5_logit <- log2((xeno_110p5+0.001)/(1-(xeno_110p5+0.001)))
testData_sub_110p5_logit <- testData_sub_110p5_logit[,order(names(testData_sub_110p5_logit))]                          

#320p2
xeno_320p2 <- testData_sub_xenograph[,c(seq(2,12,by=2),25:36)]
xeno_320p2_vpa <- xeno_320p2[,c(grep("ck",names(xeno_320p2)),grep("vpa",names(xeno_320p2)))][,c(1,2,3,7,8,9)]

testData_sub_320p2_logit <- log2((xeno_320p2+0.001)/(1-(xeno_320p2+0.001)))
testData_sub_320p2_logit <- testData_sub_320p2_logit[,order(names(testData_sub_320p2_logit))]

#test1: adaptive_B=F, adaptive_S=F, mixture_beta=F
mcmc.chain <- assign.mcmc(Y = testData_sub_TCGA_logit, Bg = B_vector, X = S_matrix, Delta_prior_p = Pi_matrix, iter=2000, adaptive_B=F, adaptive_S=F, mixture_beta=F, p_beta = 0.5)
mcmc.pos.mean1 <- assign.summary(test=mcmc.chain, burn_in=1000, iter=2000, adaptive_B=F, adaptive_S=F,mixture_beta=F)

#test2: adaptive_B=T, adaptive_S=F, mixture_beta=F
mcmc.chain <- assign.mcmc(Y = testData_sub_TCGA_logit, Bg = B_vector, X = S_matrix, Delta_prior_p = Pi_matrix, iter=2000, adaptive_B=T, adaptive_S=F, mixture_beta=F, p_beta = 0.5)
mcmc.pos.mean2 <- assign.summary(test=mcmc.chain, burn_in=1000, iter=2000, adaptive_B=T, adaptive_S=F,mixture_beta=F)

#test3: adaptive_B=T, adaptive_S=T, mixture_beta=F
mcmc.chain <- assign.mcmc(Y = testData_sub_110p5_logit, Bg = B_vector, X = S_matrix, Delta_prior_p = Pi_matrix, iter=2000, adaptive_B=T, adaptive_S=T, mixture_beta=F, p_beta = 0.5)
mcmc.pos.mean3 <- assign.summary(test=mcmc.chain, burn_in=1000, iter=2000, adaptive_B=T, adaptive_S=T,mixture_beta=F)

#test4: adaptive_B=T, adaptive_S=T, mixture_beta=T
mcmc.chain <- assign.mcmc(Y = testData_sub_110p5_logit, Bg = B_vector, X = S_matrix, Delta_prior_p = Pi_matrix, iter=2000, adaptive_B=T, adaptive_S=T, mixture_beta=T, p_beta = 0.5)
mcmc.pos.mean4 <- assign.summary(test=mcmc.chain, burn_in=1000, iter=2000, adaptive_B=T, adaptive_S=T,mixture_beta=T)
vpa_pa <- mcmc.pos.mean4$kappa_pos
row.names(vpa_pa) <- names(testData_sub_TCGA)
write.csv(vpa_pa, file="methylation_TCGA_all_500_topCor-nonadp.csv") # csv file for all TCGA tumor-normal samples


# plots and tables
label <- as.factor(c(rep("normal",32),rep("tumor",215)))
label <- as.factor(c(rep("ck",6),rep("vpa_treated",6),rep("zeb_treated",6)))

pdf("methylation_combat_TCGA_200_vpa_2h.pdf")
par(mfrow=c(2,2))
boxplot(mcmc.pos.mean1$beta ~ label,ylab="vpa signature",main="test1")
boxplot(mcmc.pos.mean2$beta ~ label,ylab="vpa signature",main="test2")
boxplot(mcmc.pos.mean3$beta ~ label,ylab="vpa signature",main="test3")
boxplot(mcmc.pos.mean4$kappa ~ label,ylab="vpa signature",main="test4")
dev.off()

pa <- matrix(mcmc.pos.mean4$kappa,nrow=ncol(testData_sub_TCGA_logit),1);rownames(pa) <- names(testData_sub_TCGA_logit)
pa1 <- pa[order(rownames(pa)),1,drop=F]
pa2 <- pa1[c(grep("11A",rownames(pa1))-1, grep("11B",rownames(pa1))-1,grep("11A",rownames(pa1)), grep("11B",rownames(pa1))),1,drop=F]
dim(pa2) <- c(32,2)
rownames(pa2) <- sapply(rownames(pa1)[1:32],function(x){substr(x,1,12)})
colnames(pa2) <- c("tumor", "normal")

write.csv(pa2,file="methylation_tcga_vpa_NEW.csv")

###
coeff <- mcmc.pos.mean4$kappa
rownames(coeff) <- names(testData_sub_logit)
colnames(coeff) <- "vpa_6h_signature"
write.csv(coeff,file="110p5_vpa6h_methylation_signature.csv")

coeff <- mcmc.pos.mean4$kappa
rownames(coeff) <- names(testData_sub_logit)
colnames(coeff) <- "vpa_6h_signature"
write.csv(coeff,file="320p2_vpa6h_methylation_signature.csv")

pdf("methylation_uni+Combat_110p5_vpa_200_6h.pdf")
boxplot(mcmc.pos.mean4$kappa ~ label,ylab="vpa signature",main="110p5_200_vpa_2h_ Unique+Combat_methylation_signature.pdf")
dev.off()                           

pdf("methylation_320p2_50_vpa_6h.pdf")
boxplot(mcmc.pos.mean4$kappa ~ label,ylab="vpa signature",main="320p2_50_vpa_6h_methylation_signature.pdf")
dev.off()                           

pdf("methylation$$_tcga_vpa6h_500_topCor.pdf")
boxplot(mcmc.pos.mean4$kappa ~ label, ylab="vpa signature",col=c("darkgreen","red"),main="TCGA_vpa_6h_methylation_topCor.pdf")
dev.off()  


## VPA methylation signature subtypes
methyl <- read.csv("VPA_Methyl_Subtypes_Supriya.csv",as.is=T)

pdf("methylation_subtype_supriya.pdf")


boxplot(methyl[,c(2,5,8,11,14,17,20)], ylab="ASSIGN Signature", las=1, at=c(1,2,4,5,6,7,8), 
        cex.lab=0.8, col="brown",frame=F, xaxt = "n", main="VPA Methylation Profile")
legend("topleft",
       legend=c("VPA"),
       col=c("brown") ,
       pch=15)
axis(1.,at=c(1,2,4,5,6,7,8),labels=c("Normal","Tumor","Basal","Her2","Lum-A","Lum-B","Normal-Like"),tick=FALSE,cex.axis=0.7)
dev.off()

############################################################


## Supriya's analysis-gene selection for gene expression by LIMMA

library(sva)

expr <- read.table("finalMerged.txt",row.names="external_gene_id",header=T)
expr2 <- expr[complete.cases(expr),]
  
### ComBat
batch <- c(rep("cellLine",36),rep("tcga",247))
expr_combat <- ComBat(dat=expr2, batch, mod=NULL,numCovs=NULL, par.prior=TRUE,prior.plots=FALSE)
write.table(expr_combat,file="exprFinal_v2_combat.txt",quote=F,sep="\t") 

cellline <- expr_combat[,1:36]
vpa.cellline <- sort(names(cellline)[c(grep("vpa",names(cellline)),grep("ctr",names(cellline)))])[-c(3,4,11,12,19)]
vpa <- cellline[,vpa.cellline]
 
vpa.matrix <- as.matrix(vpa)
group <- rep(c(1,0),length=16)
pair <- as.factor(rep(1:8,each=2,length=16))
                           

# ### select significant genes using expr_combat
result <- matrix(nrow=nrow(vpa.matrix),ncol=4)
colnames(result) <- c("est_beta0","est_beta1","tstat", "pvalue")
rownames(result) <- row.names(vpa)
for (i in 1:nrow(vpa.matrix)){
if(i%%100==0){print(i)}
lm1 <- lm(vpa.matrix[i,] ~ group + pair)
result[i,] <- c(summary(lm1)$coef[1,1],summary(lm1)$coef[2,c(1,3,4)])
}
                           
topGenes <- result[order(result[,4]),]
fdr <- p.adjust(topGenes[,4],method="fdr")
topGenes <- cbind(topGenes,fdr)
library(limma)
design <-  model.matrix(~group+pair)
fit <- lmFit(vpa.matrix,design)
fit <- eBayes(fit)                     
topGenes_vpa <-topTable(fit,coef=2,number=nrow(vpa.matrix))

                        
# 2000 gene selection
nTop <- 5000
topGenes_vpa <-topTable(fit,coef=2,number=nTop)
write.csv(topGenes_vpa, file="vpa_diffExp_GeneList_5000.csv")
                      

###########
### ASSIGN
library(ASSIGN, "/usr2/faculty/wej/R/x86_64-unknown-linux-gnu-library/2.15")
##VPA_500
geneList_vpa <- rownames(topGenes_vpa)
S_matrix <- -fit2$coefficients[geneList_vpa,2]
B_vector <- fit2$coefficients[geneList_vpa,1]+fit2$coefficients[geneList_vpa,2]
Pi_matrix <- rep(0.95,nrow(topGenes_vpa))


#TCGA
geneList <- rownames(topGenes)[1:nTop]
testData_TCGA_sub <- expr_combat[geneList,37:283]
                          
#test4: adaptive_B=T, adaptive_S=T, mixture_beta=T
mcmc.chain <- assign.mcmc(Y = testData_TCGA_sub , Bg = B_vector, X = S_matrix, Delta_prior_p = Pi_matrix, iter=2000, adaptive_B=T, adaptive_S=T, mixture_beta=T, p_beta = 0.5)
mcmc.pos.mean4 <- assign.summary(test=mcmc.chain, burn_in=1000, iter=2000, adaptive_B=T, adaptive_S=T,mixture_beta=T)
vpa_pa <- mcmc.pos.mean4$beta_pos

row.names(vpa_pa) <- names(testData_sub_TCGA)
write.csv(vpa_pa, file="expression_TCGA_all_supriya_1000.csv") # csv file for all TCGA tumor-normal samples

# plots and tables
label <- as.factor(c(rep("normal",32),rep("tumor",215)))
                           
                      
pdf("expression_TCGA_1000_supriya.pdf")
boxplot(mcmc.pos.mean4$kappa ~ label,ylab="vpa signature",main="TCGA_combat_vpa_1000_expression.pdf")
dev.off()  
                           
#############################################

## gene expression analysis with top correlated genes

library(sva)

expr <- read.table("finalMerged.txt",row.names="external_gene_id",header=T)
expr2 <- expr[complete.cases(expr),]

batch <- c(rep("cellLine",36),rep("tcga",247))
expr_combat <- ComBat(dat=expr2, batch, mod=NULL,numCovs=NULL, par.prior=TRUE,prior.plots=FALSE)

 
###########
### ASSIGN
library(ASSIGN, "/usr2/faculty/wej/R/x86_64-unknown-linux-gnu-library/2.15")
           
topGenes_unique <- read.csv("correlation_unique_500_geneList.csv", header = TRUE)

geneList <- topGenes_unique[,1]
S_matrix <- -fit2$coefficients[geneList,2]
B_vector <- fit2$coefficients[geneList,1]+fit2$coefficients[geneList,2]
Pi_matrix <- rep(0.95,length(S_matrix))                                                  
                           
#TCGA

testData_TCGA_sub <- expr_combat[geneList,37:283]
 

#test4: adaptive_B=T, adaptive_S=T, mixture_beta=T
mcmc.chain <- assign.mcmc(Y = testData_TCGA_sub , Bg = B_vector, X = S_matrix, Delta_prior_p = Pi_matrix, iter=2000, adaptive_B=F, adaptive_S=F, mixture_beta=T, p_beta = 0.5)
mcmc.pos.mean4 <- assign.summary(test=mcmc.chain, burn_in=1000, iter=2000, adaptive_B=F, adaptive_S=F,mixture_beta=T)
vpa_pa <- mcmc.pos.mean4$beta_pos
                           
row.names(vpa_pa) <- names(testData_sub_TCGA)
write.csv(vpa_pa, file="expression_TCGA_all_supriya_topCor_nonadp.csv") # csv file for all TCGA tumor-normal samples
                           
# plots and tables
label <- as.factor(c(rep("normal",32),rep("tumor",215)))
                           
pdf("expression_TCGA_1000_supriya.pdf")
boxplot(mcmc.pos.mean4$kappa ~ label,ylab="vpa signature",main="TCGA_combat_vpa_1000_expression.pdf")
dev.off() 