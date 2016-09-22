

### barplots for overlapping cell lines in ICBP with Andrea's cell line
icbp_subtype <- read.csv(file="icbp_subtype.csv",header=TRUE)
icbp<- as.matrix(icbp_subtype)
pdf("ICBP_subtype.pdf")
barplot(icbp,col=rep(c("brown","forestgreen","orange"),11),beside=TRUE,cex.axis=1,cex.names=.5,space=rep(c(1,0,0),11),
        ylim=0:1,legend=c("VPA","SAHA","TSA"),cex=0.1,
        ylab= "ASSIGN Signature")
dev.off()

### Boxplots for all the cell lines in ICBP
icbp_allCellLines <- read.csv("icbp_allCellLines.csv")
pdf("ICBP_allCellLines.pdf")
boxplot(c(icbp_allCellLines[2:4],icbp_allCellLines[6:8],icbp_allCellLines[10:12],icbp_allCellLines[14:16]),
          col=rep(c("brown","forestgreen","orange"),4),
          at=c(1,2,3,6,7,8,11,12,13,16,17,18),
        main="Multi-Drug Expression Signature Validation in ICBP Breast Cancer Cell lines", cex.main=1)
legend("topleft",c("VPA","SAHA","TSA"),fill=c("brown","forestgreen","orange"))
dev.off()

icbp_luminal <- read.csv("ICBP_luminalONLY.csv", header=TRUE)
pdf("ICBP_Luminal.pdf")
boxplot(icbp_luminal[2:4],col=c("brown","forestgreen","orange"))
dev.off()
