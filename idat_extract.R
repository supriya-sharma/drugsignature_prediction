#source("http://bioconductor.org/biocLite.R")
#biocLite("lumi")
#biocLite("methylumi")
#biocLite("IlluminaHumanMethylation450k.db")


#library(methylumIDAT)
library(methylumi)
library(lumi)
setwd("~/Desktop/projects/methylation/Illumina_Methylation450K/5790865004/")


extractDir=function(directory){
  setwd(directory)
  filenames=list.files(pattern="Grn.idat")
  cat("Processing the the following .idat files:",filenames,sep='\n')
  barcode=unlist(strsplit(filenames,"_Grn.idat"))

  mldat=methylumIDAT(barcode)

  meth=methylated(mldat)
  unmeth=unmethylated(mldat)

  n=nrow(meth)
  #n=100
  probes=rownames(meth)[1:n]

  library(IlluminaHumanMethylation450k.db)
  geneanno = IlluminaHumanMethylation450kSYMBOL
  probes=probes[probes %in% keys(geneanno)]  ### filter out probes without annotation
  meth=meth[probes,]
  unmeth=unmeth[probes,]
  beta= round( meth / (meth + unmeth + 100), 4) ## Calculate Beta values
  n=nrow(meth)

    
  ### Sequence ####
  nuanno = IlluminaHumanMethylation450kNUID
  nuids=unlist(as.list(nuanno[probes]))
  seq=sapply(nuids,id2seq)

  ###  Gene name  ####
  geneanno = IlluminaHumanMethylation450kSYMBOL
  genes=unlist(as.list(geneanno[probes]))

  ### Refseq ###
  refanno=IlluminaHumanMethylation450kREFSEQ
  refseq=unlist(as.list(refanno[probes]))

  ###  Chromsome        
  chranno = IlluminaHumanMethylation450kCHR
  chr=unlist(as.list(chranno[probes]))

  ### Distance to start of gene ####
  #distanno = IlluminaHumanMethylation450kCHRLOC
  #dist=unlist(as.list(distanno[probes]))

  ###  Location ####
  locanno=IlluminaHumanMethylation450kCPG37
  loc=unlist(as.list(locanno[probes]))

  #####
  anno=cbind(paste("chr",chr,sep=''),loc,genes,seq)

  for (i in 1:length(barcode)){
    tmp=cbind(anno[,1:2],meth[1:n,i],unmeth[1:n,i],beta[1:n,i],rownames(meth)[1:n],anno[,3:4])
    colnames(tmp)=c("chromosome","location","methylated","unmethylated","beta","probe","gene","sequence")
    tmp=tmp[order(as.character(tmp[,1]),as.numeric(tmp[,2])),]
    write.table(tmp,file=paste(barcode[i],'.rawdat',sep=''),quote=F,sep='\t',row.names=F)
  }
  cat("Success!\n")
}


dir1="~/Desktop/projects/methylation/Illumina_Methylation450K/5790865004/"
extractDir(dir1)

dir1="~/Desktop/projects/methylation/Illumina_Methylation450K/5790865027/"
extractDir(dir1)

dir1="~/Desktop/projects/methylation/Illumina_Methylation450K/5806361006/"
extractDir(dir1)

dir1="~/Desktop/projects/methylation/Illumina_Methylation450K/5806361007/"
extractDir(dir1)

dir1="~/Desktop/projects/methylation/Illumina_Methylation450K/5806361011/"
extractDir(dir1)

dir1="~/Desktop/projects/methylation/Illumina_Methylation450K/5806361048/"
extractDir(dir1)

## Get betas:

dirs = c("~/Desktop/projects/methylation/Illumina_Methylation450K/5790865004/","~/Desktop/projects/methylation/Illumina_Methylation450K/5790865027/","~/Desktop/projects/methylation/Illumina_Methylation450K/5806361006/","~/Desktop/projects/methylation/Illumina_Methylation450K/5806361007/","~/Desktop/projects/methylation/Illumina_Methylation450K/5806361011/","~/Desktop/projects/methylation/Illumina_Methylation450K/5806361048/")
chr=location=gene=probe=beta=NULL
for (j in dirs){
	print(j)
	setwd(j)
	filenames=list.files(pattern=".rawdat")
	filenames=filenames[substr(filenames,nchar(filenames)-4,nchar(filenames))!=".slam"]
	
	for (i in filenames){
		print(i)
		tmp=read.table(i,header=T)
		if (is.null(probe)){
			probe = tmp$probe
			gene = tmp$gene
			chr = tmp$chromosome
			location=tmp$location	
		}
		beta=cbind(beta, tmp$beta)
		colnames(beta)[ncol(beta)]=substr(i,1,nchar(i)-7)
	}
}

setwd("~/Desktop/projects/methylation/Illumina_Methylation450K/")
all_beta = data.frame(chromosome=chr, location, probe, gene, beta)
write.table(all_beta, file="all_betas.txt",quote=F,row.names=F) 

