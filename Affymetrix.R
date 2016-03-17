## Affymetrix.R
## Tutorial script (UCLANeurogenetics on GitHub) for analyzing an Affymetrix microarray
## GSE20295 (http://www.ncbi.nlm.nih.gov/geo/), Parkinson's Disease Multi-Region Brain study

### First, install necessary libraries and set working directory

options(stringsAsFactors = FALSE)
source("http://bioconductor.org/biocLite.R")
biocLite("affy")                              ### Use bioconductor to download new libraries

library(affy) ### Load your library every time
library(GEOquery)
library(limma)
library(WGCNA)
library(sva)
library(biomaRt)

setwd("C:/Users/Jill/Dropbox/DHGLab/Tutorials/Microarray/")

###### Steps in Analysis
### (1) Get Data
### (2) Quality Control (QC) on Raw Data 
### (3) Normalization
### (4) Batch Correction
### (5) Outlier Removal
### (6) QC on Normalized Data
### (7) Covariate Analysis
### (8) Annotate Probes
### (9) Collapse Rows
### (10) Differential Expression Analysis

###### (1) Get Data

## Two Choices For Getting Data:
## (a) Load your data directly from a source into your R session (using GEOquery R library)
## (b) Download your data to your local computer or server before loading it into your R session

### (a) Load your data directly from a source into your R session (using GEOquery R library)

## GOAL: to get datMeta and datExpr for prefrontal cortex samples only 

# getGEO() to get phenotype data for entire experiment - datMeta
gse <- getGEO("GSE20295", GSEMatrix =TRUE,getGPL=FALSE)

datMeta = pData(gse[[1]])
rownames(datMeta) = datMeta[,2]
datMeta$title = gsub(" ","_",datMeta$title)
idx = which(datMeta$source_name_ch1 == "Postmortem brain prefrontal cortex")

# getGEOSuppFiles() to get .CEL files to create affy object with expression data for entire experiment - datExpr
getGEOSuppFiles("GSE20295")
# Use software such as 7-zip to extract into "GSE20295_RAW"
# manually change GSM506039_1_1364_BA9_Pm.CEL.gz to GSM506039_1364_BA9_Pm.CEL.gz first - this was an upload error
filesPFC = paste(datMeta$geo_accession,"_",datMeta$title,".CEL.gz",sep="")[idx]
data.affy = ReadAffy(celfile.path = "./GSE20295/GSE20295_RAW", filenames = filesPFC)
datExpr = exprs(data.affy)
datMeta = datMeta[idx,]

# check ordering
GSM = rownames(pData(data.affy))
GSM = substr(GSM,1,9)
idx = match(GSM, datMeta$geo_accession)
datMeta = datMeta[idx,] 

# reformat datMeta and datExpr
datMeta = datMeta[,-c(3:7,14:36)]
datMeta$characteristics_ch1 = gsub("disease state: control","CTL",datMeta$characteristics_ch1)
datMeta$characteristics_ch1 = gsub("disease state: Parkinson's disease","PKD",datMeta$characteristics_ch1)
datMeta$characteristics_ch1.1 = gsub("gender: male","M",datMeta$characteristics_ch1.1)
datMeta$characteristics_ch1.1 = gsub("gender: female","F",datMeta$characteristics_ch1.1)
datMeta$characteristics_ch1.2 = gsub("age: ","",datMeta$characteristics_ch1.2)
datMeta$characteristics_ch1.3 = gsub("brain region: ","",datMeta$characteristics_ch1.3)
colnames(datMeta)[5:8] = c("Dx","Sex","Age","Region") 
datMeta$Dx = as.factor(datMeta$Dx)

### (b) Download your data to your local computer or server before loading it into your R session

## Although option (a) is more straightforward, sometimes with large datasets (> 500 Mb) it is better to 
## 'pre-filter', or only extract/download data that you need - in this experiment, we only want prefrontal cortex.
## Additionally, if you have your own experiment, you can simply use ReadAffy() to load your own expression data.

## You can download the prefrontal cortex expression data directly from the GEO page, 
## (http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE20295). For GSE20295_RAW.tar choose 'custom' download, and only
## select files with 'BA9' and '.CEL.gz' in the name. Download these files into a new directory 'GSE20295_PFC' 
## and load it into R from there (remember to extract the raw data first with 7-zip).

## GOAL: to get datMeta and datExpr for prefrontal cortex samples only 

# getGEO() to get phenotype data for entire experiment - datMeta
gse <- getGEO("GSE20295",GSEMatrix =TRUE,getGPL=FALSE)

datMeta = pData(gse[[1]])
rownames(datMeta) = datMeta[,2]
datMeta$title = gsub(" ","_",datMeta$title)
idx = which(datMeta$source_name_ch1 == "Postmortem brain prefrontal cortex")

# manually change GSM506039_1_1364_BA9_Pm.CEL.gz to GSM506039_1364_BA9_Pm.CEL.gz first - this was an upload error
data.affy = ReadAffy(celfile.path = "./GSE20295_PFC/GSE20295_RAW")
datExpr = exprs(data.affy)
datMeta = datMeta[idx,]

# check ordering
GSM = rownames(pData(data.affy))
GSM = substr(GSM,1,9)
idx = match(GSM, datMeta$geo_accession)
datMeta = datMeta[idx,] 

# reformat datMeta and datExpr
datMeta = datMeta[,-c(3:7,14:36)]
datMeta$characteristics_ch1 = gsub("disease state: control","CTL",datMeta$characteristics_ch1)
datMeta$characteristics_ch1 = gsub("disease state: Parkinson's disease","PKD",datMeta$characteristics_ch1)
datMeta$characteristics_ch1.1 = gsub("gender: male","M",datMeta$characteristics_ch1.1)
datMeta$characteristics_ch1.1 = gsub("gender: female","F",datMeta$characteristics_ch1.1)
datMeta$characteristics_ch1.2 = gsub("age: ","",datMeta$characteristics_ch1.2)
datMeta$characteristics_ch1.3 = gsub("brain region: ","",datMeta$characteristics_ch1.3)
colnames(datMeta)[5:8] = c("Dx","Sex","Age","Region") 
datMeta$Dx = as.factor(datMeta$Dx)

###### (2) Quality Control (QC) on Raw Data 

##primary analysis - examine distribution of expression data across samples

datExpr = log2(datExpr)
dim(datExpr)

#boxplot

boxplot(datExpr,range=0, col = as.numeric(datMeta$Dx), xaxt='n', xlab = "Array", main = "Boxplot Pre-Normalization", ylab = "Intensity")
legend("topright",legend = levels(datMeta$Dx),fill = as.numeric(as.factor(levels(datMeta$Dx))))

#histogram

i=1
plot(density((datExpr[,i]),na.rm=T),col = as.numeric(datMeta$Dx[i]),
     main = "Hist of Log2 Exp", xlab="log2 exp",xlim=c(4,16),ylim=c(0,0.5))
for(i in 2:29){
  lines(density((datExpr[,i]),na.rm=T), col = as.numeric(datMeta$Dx)[i],)
}
legend("topright",legend = levels(datMeta$Dx),fill = as.numeric(as.factor(levels(datMeta$Dx))))

#mdsplot

mds = cmdscale(dist(t(datExpr)),eig=TRUE)
plot(mds$points,col=as.numeric(datMeta$Dx),pch=19)
legend("topright",legend = levels(datMeta$Dx),fill = as.numeric(as.factor(levels(datMeta$Dx))))

###### (3) Normalization

datExpr = rma(data.affy, background=T, normalize=T, verbose=T)
datExpr = exprs(datExpr)

###### (4) Batch Correction

batch = protocolData(data.affy)$ScanDate
batch = substr(batch,1,8)
batch = as.factor(batch)
table(batch)
datMeta$Batch = batch

plot(mds$points,col=as.numeric(datMeta$Batch),pch=19)
legend("topright",legend = levels(datMeta$Batch),fill = as.numeric(as.factor(levels(datMeta$Batch))))

# remove singular batches
to_remove = (datMeta$Batch == "09/04/03")
datExpr = datExpr[,!to_remove]
datMeta = datMeta[!to_remove,]
datMeta$Batch = droplevels(datMeta$Batch)

# remove batch contributions from expression
mod = model.matrix(~datMeta$Dx)   
batch = as.factor(datMeta$Batch)
datExpr.combat = ComBat(dat = datExpr,batch = batch,mod = mod)

datExpr = datExpr.combat

mds = cmdscale(dist(t(datExpr)),eig=TRUE)
plot(mds$points,col=as.numeric(datMeta$Dx),pch=19)
legend("topright",legend = levels(datMeta$Dx),fill = as.numeric(as.factor(levels(datMeta$Dx))))

plot(mds$points,col=as.numeric(datMeta$Batch),pch=19)
legend("topright",legend = levels(datMeta$Batch),fill = as.numeric(as.factor(levels(datMeta$Batch))))

###### (5) Outlier Removal

colnames(datExpr)=datMeta$geo_accession
tree = hclust(dist(t(datExpr)), method="average")
plot(tree)

normadj = (0.5 + 0.5*bicor(datExpr))^2
netsummary = fundamentalNetworkConcepts(normadj)
C = netsummary$Connectivity
Z.C = (C-mean(C))/sqrt(var(C))
to_keep = abs(Z.C) < 2
table(to_keep)
colnames(datExpr)[!to_keep]

datExpr = datExpr[,to_keep]
datMeta = datMeta[to_keep,]

###### (6) QC on Normalized Data

dim(datExpr)
dim(datMeta)

#boxplot

boxplot(datExpr,range=0, col = as.numeric(datMeta$Dx), xaxt='n', xlab = "Array", main = "Boxplot Normalized", ylab = "Intensity")
legend("topright",legend = levels(datMeta$Dx),fill = as.numeric(as.factor(levels(datMeta$Dx))))

#histogram

i=1
plot(density((datExpr[,i]),na.rm=T),col = as.numeric(datMeta$Dx[i]),main = "Hist of Normalized Exp", xlab="log2 exp")
for(i in 2:27){
  lines(density((datExpr[,i]),na.rm=T), col = as.numeric(datMeta$Dx)[i],)
}
legend("topright",legend = levels(datMeta$Dx),fill = as.numeric(as.factor(levels(datMeta$Dx))))

#mdsplot

mds = cmdscale(dist(t(datExpr)),eig=TRUE)
plot(mds$points,col=as.numeric(datMeta$Dx),pch=19)
legend("topright",legend = levels(datMeta$Dx),fill = as.numeric(as.factor(levels(datMeta$Dx))))

###### (7) Covariate Analysis


par(mfrow=c(2,2))
par(mar=c(3,2,3,2))

plot(datMeta$Dx, ylab="Number", main="Subjects")

for(i in c(6,7,9)){
  
  if( i == 6 || i == 9 ){
    print(paste(i,"Character Graph",sep=" "))
    A = anova(lm(as.numeric(as.factor(datMeta[,i])) ~ datMeta$Dx)); p = A$"Pr(>F)"[1]   
    plot(as.factor(datMeta[,i]) ~ datMeta$Dx, main=paste(colnames(datMeta)[i]," p=", signif(p,2)), ylab="", xlab="")
  }
  else{
    print(paste(i,"Number Graph",sep=" "))
    A = anova(lm(as.numeric(datMeta[,i]) ~ datMeta$Dx)); p = A$"Pr(>F)"[1]   
    plot(as.numeric(as.character(datMeta[,i])) ~ datMeta$Dx, main=paste(colnames(datMeta)[i]," p=", signif(p,2)), ylab="", xlab="")
  }
}

###### (8) Annotate Probes

ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="www.ensembl.org")

f = listFilters(ensembl)
a = listAttributes(ensembl)

identifier <- "affy_hg_u133a"
getinfo <- c("affy_hg_u133a", "ensembl_gene_id", "entrezgene", "external_gene_name")
geneDat <- getBM(attributes = getinfo, filters=identifier, values = rownames(datExpr),mart=ensembl)

idx = match(rownames(datExpr),geneDat$affy_hg_u133a)
geneDat = geneDat[idx,]

table(is.na(geneDat$ensembl_gene_id))

to_keep = (is.na(geneDat$ensembl_gene_id) == FALSE)
geneDat = geneDat[to_keep,]
datExpr = datExpr[to_keep,]

dim(datExpr)
dim(geneDat)

###### (9) Collapse Rows

table(duplicated(geneDat$affy_hg_u133a))
table(duplicated(geneDat$ensembl_gene_id))

CR = collapseRows(datExpr, rowGroup = geneDat$ensembl_gene_id, rowID = geneDat$affy_hg_u133a)
datExpr = CR$datETcollapsed
idx = match(CR$group2row[,"selectedRowID"], geneDat$affy_hg_u133a)
geneDat = geneDat[idx,]
rownames(geneDat) = geneDat$ensembl_gene_id

dim(datExpr)
dim(geneDat)
dim(datMeta)

write.csv(datExpr, file = "DatExpr.csv")
write.csv(datMeta, file = "DatMeta.csv")
write.csv(geneDat, file = "geneDat.csv")

datExpr=read.csv("DatExpr.csv",row.names=1)
datMeta=read.csv("DatMeta.csv")
geneDat=read.csv("geneDat.csv")

###### (10) Differential Expression Analysis

datMeta$Sex=factor(datMeta$Sex,levels=c("M","F"))
mod = model.matrix(~datMeta$Dx+datMeta$Sex+as.numeric(datMeta$Age))
fit = lmFit(datExpr,mod)
fit = eBayes(fit)
tt = topTable(fit,coef = 2,n = Inf,genelist = geneDat)

par(mfrow=c(1,1))
hist(tt$P.Value)

