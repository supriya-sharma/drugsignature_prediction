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
strain <- as.factor(rep(1:6,times=3))

design <-  model.matrix(~treatment+strain)

fit2 <- lmFit(vpa_cells_logit,design)
fit2 <- eBayes(fit2)

nTop <- 5000
topGenes_vpa2h_2 <-topTable(fit2,coef=2,number=nTop)
topGenes_vpa6h_2 <- topTable(fit2, coef=3,number=nTop)

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

#TCGA

testData_sub_TCGA <-met_combat[topCor,73:319]
testData_sub_TCGA[testData_sub_TCGA<0] = 0
testData_sub_TCGA[testData_sub_TCGA>.998] = .998
testData_sub_TCGA_logit <- log2((testData_sub_TCGA+0.001)/(1-(testData_sub_TCGA+0.001)))

### ASSIGN

library(ASSIGN, "/usr2/faculty/wej/R/x86_64-unknown-linux-gnu-library/2.15")


assign.wrapper(trainingData=NULL, testData=testData_sub_TCGA_logit, trainingLabel=NULL, testLabel=NULL, geneList=topCor,n_sigGene=NULL, adaptive_B=TRUE, adaptive_S=TRUE, mixture_beta=TRUE, iter=2000, burn_in=1000)

## gene expression analysis with top correlated gene

expr <- read.table("finalMerged.txt",row.names="external_gene_id",header=T)
expr2 <- expr[complete.cases(expr),]

batch <- c(rep("cellLine",36),rep("tcga",247))
expr_combat <- ComBat(dat=expr2, batch, mod=NULL,numCovs=NULL, par.prior=TRUE,prior.plots=FALSE)

#TCGA

testData_sub_TCGA <- expr_combat[geneList,37:283]

assign.wrapper(trainingData=NULL, testData=testData_sub_TCGA, trainingLabel=NULL, testLabel=NULL, geneList=topCor, n_sigGene=NULL, adaptive_B=TRUE,
               adaptive_S=TRUE, mixture_beta=TRUE, iter=2000, burn_in=1000)

