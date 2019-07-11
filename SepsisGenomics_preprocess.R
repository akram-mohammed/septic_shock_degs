#Author : Akram Mohammed and Rishikesan Kamaleswaran
install.packages("magrittr")
install.packages("dplyr")
install.packages("AnnotatedDataFrame")
#install the core bioconductor packages, if not already installed
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite(BiocUpgrade)
# install additional bioconductor libraries, if not already installed
biocLite("GEOquery")
biocLite("affy")
biocLite("gcrma")
biocLite("hgu133plus2cdf")
biocLite("hgu133plus2probe")
biocLite("hgu133plus2.db")

library(limma)
library(magrittr)
#Load the necessary libraries
library(GEOquery)
library(affy)
library(gcrma)
library(hgu133plus2cdf)
library(hgu133plus2probe)
library(hgu133plus2.db)
library(dplyr)
library(AnnotationDbi)
#Set working directory for download
setwd("C:/Users/amoham18/Documents/SepsisGenomics")

#Download the CEL file package for this dataset (by GSE - Geo series id)
getGEOSuppFiles("GSE66099")

#Unpack the CEL files
setwd("C:/Users/amoham18/Documents/SepsisGenomics/GSE66099")

untar("GSE66099_RAW.tar", exdir="data")
cels = list.files("data/", pattern = "CEL")
sapply(paste("data", cels, sep="/"), gunzip)
cels = list.files("data/", pattern = "CEL")

setwd("C:/Users/amoham18/Documents/SepsisGenomics/GSE66099/data")

pd<-read.AnnotatedDataFrame("SepticShock_ClassLabel_SS.csv",sep=",", header=T)

rawData=ReadAffy(verbose=TRUE, filenames=cels, cdfname="hgu133plus2", phenoData=pd) #From bioconductor

#perform RMA normalization
normData=rma(rawData)

#Get the important stuff out of the data - the expression estimates for each array
rma=exprs(normData)

#Filter Affy controls
rma=rma[which(!grepl("AFFX", rownames(rma))),]

#Format values to 5 decimal places
rma=format(rma, digits=5)

#Map probe set identifiers to Entrez gene symbols and IDs and then combine with raw data.
#Map probe sets to gene symbols or other annotations
#To see all available mappings for this platform
ls("package:hgu133plus2.db") #Annotations at the exon probeset level

probes=rownames(rma)

collapser <- function(x){
  x %>% unique %>% sort %>% paste(collapse = "|")
}

# Redefinition of ae.annots
newSymbols <- AnnotationDbi::select(
  x       = hgu133plus2.db,
  keys    = rownames(rma),
  columns = c("PROBEID", "ENSEMBL", "ENTREZID", "SYMBOL"),
  keytype = "PROBEID"
) %>%
  group_by(PROBEID) %>%
  summarise_all(funs(collapser)) %>%
  ungroup

fd <- new("AnnotatedDataFrame",
          data = data.frame(rma, stringsAsFactors = FALSE)
)
rownames(fd) <- newSymbols$PROBEID

#Combine gene annotations with raw data
newrma=cbind(newSymbols,rma)

#Write RMA-normalized, mapped data to file
write.table(newrma, file = "rma2.txt", quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

# Run the local jupyternotebook to extract the input files for sepsisGenomics_DGE.R