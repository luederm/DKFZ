# @author Matthew Lueder
# @def Combines clinical and gene expression data from TCGA
# Make sure expression data was collected using HG-U133 Plus2 platform. This script uses justRMA to normalize

# Set this to base directory crreated when instaling + extracting TCGA data
setwd("C:\\Users\\Matthew\\Desktop\\DKFZ\\TCGA\\DL3")

source("http://bioconductor.org/biocLite.R") 

#biocLite("affy")
library(affy)


# Download/Install annotation packages
#biocLite("hgu133plus2cdf", suppressUpdates=T)
library(hgu133plus2cdf)   
#biocLite("hgu133plus2.db", suppressUpdates=T)
library(hgu133plus2.db)   
#biocLite("hgu133plus2probe", suppressUpdates=T)
library(hgu133plus2probe)

# Get list of CEL files
CEL_paths <- dir('./Expression-Genes/WUSM__HG-U133_Plus_2/Level_1', full.names = T, pattern = "CEL$")

# Read CEL files and Normalize/background correct 
eset <- justRMA(filenames = CEL_paths)

# Get expression data in form of a data frame
exprs <- exprs(eset) 

# Assign gene names to probe ids
sid <- rownames(exprs)
sym <- unlist(mget(sid, hgu133plus2SYMBOL, ifnotfound = NA))
rownames(exprs) <- sym
rm(sid, sym)

# Create patient barcode feild from CEL file name
toPatientBarcodes <- function(CEL.names) {
  nameVec <- NULL
  for (name in CEL.names) {
    index <- gregexpr(patter = '-', name)[[1]][3] - 1
    nameVec <- c(nameVec, substr(name, 0, index))
  }
  return(nameVec)
}
# Changes precision: How to prevent?
combined <- rbind( exprs, "bcr_patient_barcode" = toPatientBarcodes(colnames(exprs)) )

# Read in clinical data
clinical <- read.table("./Clinical/Biotab/nationwidechildrens.org_clinical_patient_laml.txt",
                sep = "\t", header = T, stringsAsFactors = F)

# Reduce to patient information TODO: Keep CDE identifiers in separate object first?
clinical <- clinical[3:nrow(clinical),]

# Recode values
clinical[clinical == "NO"] <- 0
clinical[clinical == "No"] <- 0

clinical[clinical == "YES"] <- 1
clinical[clinical == "Yes"] <- 1

clinical[clinical == "[Not Available]"] <- NA
clinical[clinical == "null"] <- NA

# Set up status and time feilds for survival analysis
status <- NULL
time <- NULL
for (i in 1:nrow(clinical)) {
  if (clinical[i,"vital_status"] == "Dead") {
    status <- c(status, 1)
    time <- c(time, clinical[i,"death_days_to"])
  }
  else {
    status <- c(status, 0)
    time <- c(time, clinical[i,"last_contact_days_to"])
  }
}
clinical$status <- status
clinical$time <- time

combined <- merge( t(combined), clinical )

colnames(clinical)














