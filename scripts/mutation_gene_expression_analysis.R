require(monocle)
#library(ggplot2)
#library(reshape2)
require(SC3)
require(argpars)
require(scater)
require(readr)
require(stringr)
suppressPackageStartupMessages(library("argparse"))

parser <- ArgumentParser()

# specify our desired options 
# by default ArgumentParser will add an help option 
parser$add_argument("-s", "--sample-sheet", default="CuffNorm/samples.table",
                    help="CuffNorm sample sheet")
parser$add_argument("-g", "--gene-annotation", default="CuffNorm/genes.attr_table",
                    help="Gene annotation used for CuffNorm")
parser$add_argument("-f", "--fpkm-matrix", default="CuffNorm/genes.fpkm_table", 
                    help="CuffNorm FPKM matrix")
parser$add_argument("-a", "--annotation", default="annotation.csv", 
                    help="Annotation file for each cell. First column contains cell ID")

args <- parser$parse_args()
#print("reading input files...")
#print(args)
#read CuffNorm output
sample_sheet <- read.delim(args$sample_sheet, row.names=1)
gene_annotation <- read.delim(args$gene_annotation, row.names=1)
fpkm_matrix <- read.delim(args$fpkm_matrix, row.names=1)


#make sure fpkm and sample sheet has the same Cell ID
colnames(fpkm_matrix)<-rownames(sample_sheet)

#read annotation
annotation <- read.csv(args$annotation, header=TRUE, row.names=1) 
rownames(annotation)<-paste(rownames(annotation),"__0",sep="") #add __0 that is added by VDJPuzzle

#merge annotation to sample sheet
sample_sheet<-merge(sample_sheet,annotation, by="row.names", all.x = T)
rownames(sample_sheet)<-sample_sheet$Row.names
sample_sheet<-sample_sheet[colnames(fpkm_matrix),]

#create CellDataSet object
pd <- new("AnnotatedDataFrame", data = sample_sheet)
fd <- new("AnnotatedDataFrame", data = gene_annotation)
my_data <- newCellDataSet(as.matrix(fpkm_matrix), phenoData = pd, featureData = fd)


#filter lowly expressed genes
min_exp10 <- detectGenes(my_data, min_exp = 10)
expressed_genes10 <- subset(fData(min_exp10), num_cells_expressed >= 1)
min_exp10 <- min_exp10[rownames(expressed_genes10),]

#remove ERCC if any
min_exp10 <- min_exp10[!rownames(min_exp10) %in% rownames(min_exp10)[grep("ERCC-", rownames(min_exp10))], ]


### sc3 and mutations ###
#read mutations from VDJPuzzle output and normalize by region length
IGH.vdjpuzzle <- read_delim("summary_corrected/igh.csv", "\t", escape_double = FALSE, trim_ws = TRUE)
IGK.vdjpuzzle <- read_delim("summary_corrected/igk.csv", "\t", escape_double = FALSE, trim_ws = TRUE)
IGL.vdjpuzzle <- read_delim("summary_corrected/igl.csv", "\t", escape_double = FALSE, trim_ws = TRUE)

IGH.vdjpuzzle$CellID<-unlist(lapply(IGH.vdjpuzzle$CellID,function(x) strsplit(x,"IGH_")[[1]][2]))
IGK.vdjpuzzle$CellID<-unlist(lapply(IGK.vdjpuzzle$CellID,function(x) strsplit(x,"IGK_")[[1]][2]))
IGL.vdjpuzzle$CellID<-unlist(lapply(IGL.vdjpuzzle$CellID,function(x) strsplit(x,"IGL_")[[1]][2]))
mutations.h<-c()
avgPosCDR3<-c()
avgPosCDR2<-c()
for (i in 1:dim(IGH.vdjpuzzle)[1]) {
  mfr1<-IGH.vdjpuzzle[i,"mutations.nt.FR1"]
  mfr2<-IGH.vdjpuzzle[i,"mutations.nt.FR2"]
  mfr3<-IGH.vdjpuzzle[i,"mutations.nt.FR3"]
  mfr4<-IGH.vdjpuzzle[i,"mutations.nt.FR4"]
  mcdr1<-IGH.vdjpuzzle[i,"mutations.nt.CDR1"]
  mcdr2<-IGH.vdjpuzzle[i,"mutations.nt.CDR2"]
  mcdr3<-IGH.vdjpuzzle[i,"mutations.nt.CDR3"]
  
  mfr1<-sum(!is.na(strsplit(as.character(mfr1),split = ",")[[1]]))/nchar(IGH.vdjpuzzle[i,"fr1nt"])
  mfr2<-sum(!is.na(strsplit(as.character(mfr2),split = ",")[[1]]))/nchar(IGH.vdjpuzzle[i,"fr2nt"])
  mfr3<-sum(!is.na(strsplit(as.character(mfr3),split = ",")[[1]]))/nchar(IGH.vdjpuzzle[i,"fr3nt"])
  mfr4<-sum(!is.na(strsplit(as.character(mfr4),split = ",")[[1]]))/nchar(IGH.vdjpuzzle[i,"fr4nt"])
  mcdr1<-sum(!is.na(strsplit(as.character(mcdr1),split = ",")[[1]]))/nchar(IGH.vdjpuzzle[i,"cdr1nt"])
  mcdr2<-sum(!is.na(strsplit(as.character(mcdr2),split = ",")[[1]]))/nchar(IGH.vdjpuzzle[i,"cdr2nt"])
  mcdr3<-sum(!is.na(strsplit(as.character(mcdr3),split = ",")[[1]]))/nchar(IGH.vdjpuzzle[i,"cdr3nt"])
  
  mut<-c(mfr1,mfr2,mfr3,mfr4,mcdr1,mcdr2,mcdr3)
  mutations.h<-rbind(mutations.h,mut)
  
  if (mcdr3>0) {
    tmp<-strsplit(as.character(IGH.vdjpuzzle[i,"mutations.nt.CDR3"]),split = ",")[[1]]
    fromVstart<-unlist(lapply(tmp,function(x) as.numeric(substring(str_split(x,":")[[1]][1],2))))
    fromCDR3start<-(fromVstart-(nchar(IGH.vdjpuzzle[i,"fr1nt"])+nchar(IGH.vdjpuzzle[i,"fr2nt"])+nchar(IGH.vdjpuzzle[i,"fr3nt"])+nchar(IGH.vdjpuzzle[i,"cdr1nt"])+nchar(IGH.vdjpuzzle[i,"cdr2nt"])))/nchar(IGH.vdjpuzzle[i,"cdr3nt"])
    avgPosCDR3<-c(avgPosCDR3,fromCDR3start)
  }
  
  
  if (mcdr2>0) {
    tmp<-strsplit(as.character(IGH.vdjpuzzle[i,"mutations.nt.CDR2"]),split = ",")[[1]]
    fromVstart<-unlist(lapply(tmp,function(x) as.numeric(substring(str_split(x,":")[[1]][1],2))))
    fromCDR2start<-(fromVstart-(nchar(IGH.vdjpuzzle[i,"fr1nt"])+nchar(IGH.vdjpuzzle[i,"fr2nt"])+nchar(IGH.vdjpuzzle[i,"cdr1nt"])))/nchar(IGH.vdjpuzzle[i,"cdr2nt"])
    avgPosCDR2<-c(avgPosCDR2,fromCDR2start)
  }
}
rownames(mutations.h)<-IGH.vdjpuzzle$CellID
colnames(mutations.h)<-c("fr1","fr2","fr3","fr4","cdr1","cdr2","cdr3")


avgPosCDR3<-c()
mutations.l<-c()
for (i in 1:dim(IGL.vdjpuzzle)[1]) {
  mfr1<-IGL.vdjpuzzle[i,"mutations.nt.FR1"]
  mfr2<-IGL.vdjpuzzle[i,"mutations.nt.FR2"]
  mfr3<-IGL.vdjpuzzle[i,"mutations.nt.FR3"]
  mfr4<-IGL.vdjpuzzle[i,"mutations.nt.FR4"]
  mcdr1<-IGL.vdjpuzzle[i,"mutations.nt.CDR1"]
  mcdr2<-IGL.vdjpuzzle[i,"mutations.nt.CDR2"]
  mcdr3<-IGL.vdjpuzzle[i,"mutations.nt.CDR3"]
  
  mfr1<-sum(!is.na(strsplit(as.character(mfr1),split = ",")[[1]]))/nchar(IGL.vdjpuzzle[i,"fr1nt"])
  mfr2<-sum(!is.na(strsplit(as.character(mfr2),split = ",")[[1]]))/nchar(IGL.vdjpuzzle[i,"fr2nt"])
  mfr3<-sum(!is.na(strsplit(as.character(mfr3),split = ",")[[1]]))/nchar(IGL.vdjpuzzle[i,"fr3nt"])
  mfr4<-sum(!is.na(strsplit(as.character(mfr4),split = ",")[[1]]))/nchar(IGL.vdjpuzzle[i,"fr4nt"])
  mcdr1<-sum(!is.na(strsplit(as.character(mcdr1),split = ",")[[1]]))/nchar(IGL.vdjpuzzle[i,"cdr1nt"])
  mcdr2<-sum(!is.na(strsplit(as.character(mcdr2),split = ",")[[1]]))/nchar(IGL.vdjpuzzle[i,"cdr2nt"])
  mcdr3<-sum(!is.na(strsplit(as.character(mcdr3),split = ",")[[1]]))/nchar(IGL.vdjpuzzle[i,"cdr3nt"])
  
  mut<-c(mfr1,mfr2,mfr3,mfr4,mcdr1,mcdr2,mcdr3)
  mutations.l<-rbind(mutations.l,mut)
  
  if (mcdr3>0) {
    tmp<-strsplit(as.character(IGL.vdjpuzzle[i,"mutations.nt.CDR3"]),split = ",")[[1]]
    fromVstart<-unlist(lapply(tmp,function(x) as.numeric(substring(str_split(x,":")[[1]][1],2))))
    fromCDR3start<-(fromVstart-(nchar(IGL.vdjpuzzle[i,"fr1nt"])+nchar(IGL.vdjpuzzle[i,"fr2nt"])+nchar(IGL.vdjpuzzle[i,"fr3nt"])+nchar(IGL.vdjpuzzle[i,"cdr1nt"])+nchar(IGL.vdjpuzzle[i,"cdr2nt"])))/nchar(IGL.vdjpuzzle[i,"cdr3nt"])
    avgPosCDR3<-c(avgPosCDR3,fromCDR3start)
  }
  
}
rownames(mutations.l)<-IGL.vdjpuzzle$CellID
colnames(mutations.l)<-c("fr1","fr2","fr3","fr4","cdr1","cdr2","cdr3")


#avgPosCDR3<-c()
mutations.k<-c()
for (i in 1:dim(IGK.vdjpuzzle)[1]) {
  mfr1<-IGK.vdjpuzzle[i,"mutations.nt.FR1"]
  mfr2<-IGK.vdjpuzzle[i,"mutations.nt.FR2"]
  mfr3<-IGK.vdjpuzzle[i,"mutations.nt.FR3"]
  mfr4<-IGK.vdjpuzzle[i,"mutations.nt.FR4"]
  mcdr1<-IGK.vdjpuzzle[i,"mutations.nt.CDR1"]
  mcdr2<-IGK.vdjpuzzle[i,"mutations.nt.CDR2"]
  mcdr3<-IGK.vdjpuzzle[i,"mutations.nt.CDR3"]
  
  mfr1<-sum(!is.na(strsplit(as.character(mfr1),split = ",")[[1]]))/nchar(IGK.vdjpuzzle[i,"fr1nt"])
  mfr2<-sum(!is.na(strsplit(as.character(mfr2),split = ",")[[1]]))/nchar(IGK.vdjpuzzle[i,"fr2nt"])
  mfr3<-sum(!is.na(strsplit(as.character(mfr3),split = ",")[[1]]))/nchar(IGK.vdjpuzzle[i,"fr3nt"])
  mfr4<-sum(!is.na(strsplit(as.character(mfr4),split = ",")[[1]]))/nchar(IGK.vdjpuzzle[i,"fr4nt"])
  mcdr1<-sum(!is.na(strsplit(as.character(mcdr1),split = ",")[[1]]))/nchar(IGK.vdjpuzzle[i,"cdr1nt"])
  mcdr2<-sum(!is.na(strsplit(as.character(mcdr2),split = ",")[[1]]))/nchar(IGK.vdjpuzzle[i,"cdr2nt"])
  mcdr3<-sum(!is.na(strsplit(as.character(mcdr3),split = ",")[[1]]))/nchar(IGK.vdjpuzzle[i,"cdr3nt"])
  
  mut<-c(mfr1,mfr2,mfr3,mfr4,mcdr1,mcdr2,mcdr3)
  mutations.k<-rbind(mutations.k,mut)
  
  
  if (mcdr3>0) {
    tmp<-strsplit(as.character(IGK.vdjpuzzle[i,"mutations.nt.CDR3"]),split = ",")[[1]]
    fromVstart<-unlist(lapply(tmp,function(x) as.numeric(substring(str_split(x,":")[[1]][1],2))))
    fromCDR3start<-(fromVstart-(nchar(IGK.vdjpuzzle[i,"fr1nt"])+nchar(IGK.vdjpuzzle[i,"fr2nt"])+nchar(IGK.vdjpuzzle[i,"fr3nt"])+nchar(IGK.vdjpuzzle[i,"cdr1nt"])+nchar(IGK.vdjpuzzle[i,"cdr2nt"])))/nchar(IGK.vdjpuzzle[i,"cdr3nt"])
    avgPosCDR3<-c(avgPosCDR3,fromCDR3start)
  }
}
rownames(mutations.k)<-IGK.vdjpuzzle$CellID
colnames(mutations.k)<-c("fr1","fr2","fr3","fr4","cdr1","cdr2","cdr3")


#merge mutation with phenotype data
#pData(min_exp10)$CellID<-unlist(lapply(pData(min_exp10)$Row.names, function(x) strsplit(x,split = "-")[[1]][1])) # I need to remove this

mutations.h.new<-c()
for (c in unique(rownames(mutations.h))) {
  if (is.numeric(dim(mutations.h[rownames(mutations.h)==c,])[1])) {
    mutations.h.new<-rbind(mutations.h.new,colMeans(mutations.h[rownames(mutations.h)==c,]))
  }
  else {
    mutations.h.new<-rbind(mutations.h.new,mutations.h[rownames(mutations.h)==c,])
  }
}
rownames(mutations.h.new)<-unique(rownames(mutations.h))

mutations.light<-rbind(mutations.l,mutations.k)
mutations.l.new<-c()
for (c in unique(rownames(mutations.light))) {
  if (is.numeric(dim(mutations.light[rownames(mutations.light)==c,])[1])) {
    mutations.l.new<-rbind(mutations.l.new,colMeans(mutations.light[rownames(mutations.light)==c,]))
  }
  else {
    mutations.l.new<-rbind(mutations.l.new,mutations.light[rownames(mutations.light)==c,])
  }
}
rownames(mutations.l.new)<-unique(rownames(mutations.light))

class(mutations.h.new) <- "numeric" 
class(mutations.l.new) <- "numeric" 
mutations.h.new<-data.frame(mutations.h.new)
mutations.l.new<-data.frame(mutations.l.new)

colnames(mutations.h.new)<-unlist(lapply(colnames(mutations.h.new), function(x) paste0("IGH.",x)))
colnames(mutations.l.new)<-unlist(lapply(colnames(mutations.l.new), function(x) paste0("IGL.",x)))
mutations.new<-merge(mutations.h.new,mutations.l.new,by="row.names",all=TRUE)
mutations.new$Row.names<-paste(mutations.new$Row.names,"_0",sep="") #add __0 that is added by VDJPuzzle
rownames(mutations.new)<-mutations.new$Row.names

pData(min_exp10)<-merge(pData(min_exp10),mutations.new, by = "row.names")
#rownames(pData(min_exp10))<-pData(min_exp10)$Row.names
pData(min_exp10)[is.na(pData(min_exp10))]<-0

#run SC3
ann<-pData(min_exp10)
rownames(min_exp10)<-make.names(c(as.character(fData(min_exp10)[,"gene_short_name"])), unique=TRUE)
pd <- new("AnnotatedDataFrame", data = ann)
fd <- new("AnnotatedDataFrame", data = fData(min_exp10))
tmp <- exprs(min_exp10)
colnames(tmp) <- rownames(ann)
sceset <- newSCESet(fpkmData = tmp, phenoData = pd, featureData = fd, logExprsOffset = 1)

is_exprs(sceset) <- exprs(sceset) > 0
sceset <- calculateQCMetrics(sceset)

sceset <- sc3(sceset, ks = 2:5, biology = TRUE, n_cores = 4)

pData(sceset)$IGL.combined<-rowMeans(pData(sceset)[,c("IGL.fr1","IGL.fr2","IGL.fr3","IGL.cdr1","IGL.cdr2","IGL.cdr3")])
pData(sceset)$IGH.combined<-rowMeans(pData(sceset)[,c("IGH.fr1","IGH.fr2","IGH.fr3","IGH.cdr1","IGH.cdr2","IGH.cdr3")])

#sc3_interactive(sceset)
pdf("mutation_gene-expression_heatmap.pdf",width=10,height=8,paper='special')
sc3_plot_markers(sceset, k = 4, show_pdata=c(colnames(annotation),"IGH.combined","IGL.combined"))
dev.off()


