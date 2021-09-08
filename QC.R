library("BiocParallel")
library("parallel")
library("NOISeq")

#Grupos
factorsMM <- rnas.raw$grupo
factorsMM <- as.data.frame(factorsMM)
factorsMM <- cbind(factorsMM, rnas.raw$sample)
rownames(factorsMM) <- rnas.raw$sample
colnames(factorsMM) <- c("grupo", "samples")

#inter <- intersect(rownames(rawAML), annot$ensembl_gene_id)
#rnas.bf <- rawAML[rownames(rawAML) %in% inter,]
#annot.bf <- annot[annot$ensembl_gene_id %in% inter,]

#row.names(rnas.bf)<-rowData(rnas.bf)$ensembl_gene_id
row.names(annot)<-annot$Gene.name
#row.names(factorsB)<-factorsB$samples

mydata <- NOISeq::readData(
  data = rnas.after.3,
  factors = factorsMM, 
  length = annot[,c("Gene.name", "Length")],
  biotype = annot[,c("Gene.name", "Type")], 
  chromosome = annot[,c("Chr", "Start", "End")], 
  gc = annot[, c("Gene.name", "GC")])

w <- 1024
h <- 1024
p <- 24

# Biodetection plot. Per group.
mybiodetection <- dat(mydata, type="biodetection", factor = "grupo", k=0)
png(filename="AFTER_biodetection.Rd_%03d.png",  width=w, height=h, pointsize=p)
explo.plot(mybiodetection)

mybiodetection <- dat(mydata, type="biodetection", factor = "grupo", k=0)
png(filename="2AFTER_biodetection.Rd_%03d.png",  width=w, height=h, pointsize=p)
explo.plot(mybiodetection, samples = c(1, 2), toplot = "protein_coding", plottype = "comparison")
# Count distribution per biotype. Using count per million, only for one sample

mycountsbio <- dat(mydata, factor = NULL, type = "countsbio")
png(filename="AFTER_countsbio.png",  width=w, height=h, pointsize=p)
explo.plot(mycountsbio, toplot = 1, samples = 1, plottype = "boxplot")

## Count distribution per sample
mycountsbio <- dat(mydata, factor = NULL, type = "countsbio")
png(filename =  "AFTER_protein_coding_boxplot.png", width=w*2, height=h, pointsize=p)
explo.plot(mycountsbio, toplot = "protein_coding",
           samples = NULL, plottype = "boxplot")

png(filename =  "AFTER_protein_coding_barplot.png", width=w*2, height=h, pointsize=p)
explo.plot(mycountsbio, toplot = "protein_coding",
           samples = NULL, plottype = "barplot")

mycountsbio <- dat(mydata, factor = "grupo", type = "countsbio")

## Count distribution per Experimental factors
png(filename= "AFTER_protein_coding_boxplot_group.png", width=w, height=h, pointsize=p)
explo.plot(mycountsbio, toplot = "protein_coding",
           samples = NULL, plottype = "boxplot")

png(filename= "AFTER_protein_coding_barplot_group.png", width=w, height=h, pointsize=p)
explo.plot(mycountsbio, toplot = "protein_coding",
           samples = NULL, plottype = "barplot")

## Saturation plot. 

mysaturation <- dat(mydata, k = 0, ndepth = 7, type = "saturation")
png(filename= "AFTER_saturation.png", width=w, height=h, pointsize=p)
explo.plot(mysaturation, toplot="protein_coding",
           samples = c(1,3), yleftlim = NULL, yrightlim = NULL)

## Length bias detection factor = "grupos"
mylengthbias <- dat(mydata, type="lengthbias", factor = "grupo")
png(filename= "AFTER_Lengthbias.png", width=w, height=h, pointsize=p)
explo.plot(mylengthbias, toplot = "global")

##GC bias factor = "grupos"
mygcbias <- dat(mydata,  type="GCbias", factor = "grupo")
png(filename="AFTER_GCbias.png",  width=w, height=h, pointsize=p)
explo.plot(mygcbias, toplot = "global")

## RNA composition
mycomp <- dat(mydata, type="cd")
png(filename= "AFTER_RNAComposition.png", width=w, height=h, pointsize=p)
explo.plot(mycomp, samples=1:12)
