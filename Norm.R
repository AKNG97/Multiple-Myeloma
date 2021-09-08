library(SummarizedExperiment)
library(TCGAbiolinks)
require(EDASeq)
require(dplyr)
require(NOISeq)
library(DESeq2)

qry.rna <- GDCquery(project = "MMRF-COMMPASS",
                    data.category= "Transcriptome Profiling",
                    data.type = "Gene Expression Quantification",
                    workflow.type = "HTSeq - Counts",)

GDCdownload(qry.rna)

dat <- qry.rna[[1]][[1]]
table(as.factor(dat$sample_type))

rnas.raw <- GDCprepare(qry.rna, summarizedExperiment = TRUE)

data <- assay(rnas.raw)
rownames(data) <- rowData(rnas.raw)$external_gene_name
rownames(rnas.raw) <- rowData(rnas.raw)$external_gene_name
head(rownames(data))

# Columna Grupo
bool <- rnas.raw$sample_type == "Primary Blood Derived Cancer - Bone Marrow"
rnas.raw$grupo[bool] <- "MM-BM"
bool <- rnas.raw$sample_type == "Recurrent Blood Derived Cancer - Bone Marrow"
rnas.raw$grupo[bool] <- "MM-rBM"
rnas.raw$grupo <- as.factor(rnas.raw$grupo)
rnas.raw <- rnas.raw[,!is.na(rnas.raw$grupo)]

dim(data)
dataFilt <- TCGAanalyze_Filtering(tabDF = data,
                                  method = "quantile",
                                  qnt.cut = 0.25)
threshold <- round(dim(data)[2]/2)
ridx <- rowSums(dataFilt == 0) <= threshold
dataFilt <- dataFilt[ridx, ]
print(dim(dataFilt))    
ridx <- rowMeans(dataFilt) >= 10
dataFilt <- dataFilt[ridx, ]
print(dim(dataFilt))
rnas.raw <- rnas.raw[rownames(rnas.raw) %in% rownames(dataFilt), ]
rnas.filt <- rnas.raw[!duplicated(rownames(rnas.raw)),]

annot.raw<-read.delim(file="mart_export.txt", sep="\t")
names(annot.raw)<-c("Gene.name", "Chr", "Start", "End", "GC", "Type", "ensembl_gene_id")
annot.raw$Length <- abs(annot.raw$End - annot.raw$Start)
inter <- intersect(rownames(rnas.filt), annot.raw$Gene.name)
rnas.before <- rnas.filt[rownames(rnas.filt) %in% inter,]
annot.raw <- annot.raw[annot.raw$Gene.name %in% inter,]
annot <- annot.raw[!duplicated(annot.raw$Gene.name),]

#Norm 1
ln.data <- withinLaneNormalization(assay(rnas.before), annot$Length, which = "full")
gcn.data <- withinLaneNormalization(ln.data , annot$GC, which = "full")
norm.counts <- tmm(gcn.data, long = 1000, lc = 0, k = 0)
noiseqData <- NOISeq::readData(norm.counts, factors = as.data.frame(rnas.raw$grupo))
mydata2corr1 = NOISeq::ARSyNseq(noiseqData, norm = "n",  logtransf = TRUE)
rnas.after <- exprs(noiseqData)

#Norm 2
ln.data <- withinLaneNormalization(assay(rnas.before), annot$Length, which = "full")
gcn.data <- withinLaneNormalization(ln.data , annot$GC, which = "full")
Btwn.Norm <- betweenLaneNormalization(gcn.data, which = "full") 
norm.counts <- tmm(Btwn.Norm, long = 1000, lc = 0, k = 0)
noiseqData <- NOISeq::readData(norm.counts, factors = as.data.frame(rnas.raw$grupo))
mydata2corr1 = NOISeq::ARSyNseq(noiseqData, norm = "n",  logtransf = TRUE)
rnas.after.2 <- exprs(noiseqData)

#Norm 3
ln.data <- withinLaneNormalization(assay(rnas.before), annot$Length, which = "upper")
gcn.data <- withinLaneNormalization(ln.data , annot$GC, which = "upper")
Btwn.Norm <- betweenLaneNormalization(gcn.data, which = "full") 
norm.counts <- tmm(Btwn.Norm, long = 1000, lc = 0, k = 0)
noiseqData <- NOISeq::readData(norm.counts, factors = as.data.frame(rnas.raw$grupo))
mydata2corr1 = NOISeq::ARSyNseq(noiseqData, norm = "n",  logtransf = TRUE)
rnas.after.3 <- exprs(noiseqData)

#Norm 4
ln.data <- withinLaneNormalization(assay(rnas.before), annot$Length, which = "loess")
gcn.data <- withinLaneNormalization(ln.data , annot$GC, which = "loess")
Btwn.Norm <- betweenLaneNormalization(gcn.data, which = "full") 
norm.counts <- tmm(Btwn.Norm, long = 1000, lc = 0, k = 0)
noiseqData <- NOISeq::readData(norm.counts, factors = as.data.frame(rnas.raw$grupo))
mydata2corr1 = NOISeq::ARSyNseq(noiseqData, norm = "n",  logtransf = TRUE)
rnas.after.4 <- exprs(noiseqData)

#RNAs by Sample Type
MM_BM <- rnas1[, rnas1$sample_type == "Primary Blood Derived Cancer - Bone Marrow"]
MM_rBM <- rnas1[, rnas1$sample_type == "Recurrent Blood Derived Cancer - Bone Marrow"]

saveRDS(MM_BM, file="MM_BM.RDS")
saveRDS(MM_rBM, file="MM_BM.RDS")
saveRDS(rnas.before, file="MM_rnasBefore.RDS")
write.table(annot,"MM_annot.txt",sep="\t")
MM_annot <- read.delim("MM_annot.txt")

#PCAs
mydata.after = readData(rnas.after, factors = as.data.frame(rnas.raw$grupo))
myPCA.after = dat(mydata.after, type = "PCA")
explo.plot(myPCA.after)

mydata.after2 = readData(rnas.after.2, factors = as.data.frame(rnas.raw$grupo))
myPCA.after2 = dat(mydata.after2, type = "PCA")
explo.plot(myPCA.after2)

mydata.after3 = readData(rnas.after.3, factors = as.data.frame(rnas.raw$grupo))
myPCA.after3 = dat(mydata.after3, type = "PCA")
explo.plot(myPCA.after3)

mydata.after4 = readData(rnas.after.4, factors = as.data.frame(rnas.raw$grupo))
myPCA.after4 = dat(mydata.after4, type = "PCA")
explo.plot(myPCA.after4)

bool <- rnas$sample_type == "Primary Blood Derived Cancer - Bone Marrow"
rnas$grupo[bool] <- "MM-BM"
bool <- rnas$sample_type == "Recurrent Blood Derived Cancer - Bone Marrow"
rnas$grupo[bool] <- "MM-rBM"
rnas$grupo <- as.factor(rnas$grupo)
rnas <- rnas[,!is.na(rnas$grupo)]
dim(assay(rnas))

