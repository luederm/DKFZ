# This script performs survival analysis as was done in the NEJM study

setwd("C:\\Users\\Matthew\\Desktop\\DKFZ\\TCGA\\NEJM_paper")

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
CEL_paths <- dir('./TCGA_LAML_EXP/genome.wustl.edu_LAML.HG-U133_Plus_2.Level_1.1.4.0/', full.names = T, pattern = "CEL$")

# Read CEL files and Normalize/background correct 
eset <- justRMA(filenames = CEL_paths)

# Get expression data in form of a data frame
exprs <- exprs(eset) 

# Assign gene names to probe ids
sid <- rownames(exprs)
sym <- unlist(mget(sid, hgu133plus2SYMBOL, ifnotfound = NA))
rownames(exprs) <- sym
rm(sid, sym)

# Extract patient identifiers from file name
toPatientIdent <- function(CEL.names) {
  nameVec <- NULL
  for (name in CEL.names) {
    greg <- gregexpr(patter = '-', name)[[1]]
    start <- greg[2] + 1
    end <- greg[3] - 1
    nameVec <- c(nameVec, substr(name, start, end))
  }
  return(nameVec)
}
combined <- as.data.frame(t(exprs))
combined$TCGA.Patient.ID <- toPatientIdent(colnames(exprs)) 

#install.packages("openxlsx")
require(openxlsx)

# Read in supplementary table 1
supTbl1 <- read.xlsx("SuppTable01.xlsx", 1)
#supTbl1 <- supTbl1[1:200,1:30]
combined <- supTbl1[1:200,]

# Merge expression and survival data
#combined <- merge(supTbl1[,1:30], combined, by = "TCGA.Patient.ID", all = F)   Next line of code has same effect, but is much faster
#combined <- cbind( supTbl1[match(combined$TCGA.Patient.ID, supTbl1$TCGA.Patient.ID), ][setdiff(colnames(supTbl1), colnames(combined))], combined )

# Next line of code will remove patients with undetermined risk
#combined <- combined[combined$`RISK.(Cyto)` != "N.D.",]

# Reformat for survival analysis
model.data <- NULL

death <- combined[,"Expired?.3.31.12"] == "*"
DCFirst <- combined[,"OS.months.3.31.12"] == combined[,"EFS.months.3.31.12"] # Does death/censoring occur before relapse
deathFirst <- death & DCFirst # Died before chace of event

relapse <- !is.na(as.numeric(combined[,"Days.from.collection.to.first.Relapse"]))
# NOTE: Relapse and death occured at same time with one patient

# 0 = right censored, 1 = relapse, 2 = death before relapse
model.data$status = as.numeric(relapse)
model.data$status[deathFirst] <- 1 


model.data$time <- combined[,"EFS.months.3.31.12"]
model.data$age.over.60 <- as.numeric(combined$Age > 60)
model.data$wbc.over.16 <- as.numeric(combined$WBC > 16) #TODO: LOG Transform
model.data$cyto.good <- as.numeric(combined$`RISK.(Cyto)` == "Good")
model.data$cyto.poor <- as.numeric(combined$`RISK.(Cyto)` == "Poor")
model.data$cyto <- factor(combined$`RISK.(Cyto)`, c("Intermediate", "Good", "Poor") )

# For When using gene expression rather than mutations
# model.data$TP53 <- combined$TP53
# model.data$FLT3 <- combined$FLT3
# model.data$DNMT3A <- combined$DNMT3A

# When combined == supTbl1, This is a binary descibing if there is a mutation
model.data$TP53 <- as.numeric( !is.na(combined$TP53) )
model.data$DNMT3A <- as.numeric( !is.na(combined$DNMT3A) )
model.data$FLT3 <- as.numeric( !is.na(combined$FLT3) )
model.data$PML.RARA <- as.numeric( combined$`PML-RARA` != " " )
model.data$MYH11.CBFB <- as.numeric( combined$`MYH11-CBFB` != " " )
model.data$RUNX1.RUNX1T1 <- as.numeric( combined$`RUNX1-RUNX1T1` != " " )
model.data$NUP98.NSD1 <- as.numeric( combined$`NUP98-NSD1` != " " )

# For Overall survival model
model.data$status2 <- as.numeric(death)
model.data$time2 <- combined[,"OS.months.3.31.12"]

# For Relapse Free Survival
model.data$status3 <- model.data$status
model.data$status3[deathFirst] <- 2



# ** Event Free Survival **
require("survival")
# KM
efsKM <- survfit(Surv(time, status) ~ cyto, type="kaplan-meier", data = model.data)
plot(efsKM, main="Event-free Survival by Cytogenetic Risk Group", xlab="Event-free Survival (Months)", 
     ylab="Survival Probability", col = c(2, 4, 3), lty = c(5, 1, 4), lwd = 2)
legend("topright", legend = c("Intermed", "Good", "Poor"), lty = c(5, 1, 4), col = c(2, 4, 3), lwd = 2)

# Basic Model
coxph(Surv(time, status) ~ cyto + wbc.over.16 + strata(age.over.60), data = model.data)

# Add gene terms separately
geneMuts <- c("TP53", "DNMT3A", "FLT3", "PML.RARA", "MYH11.CBFB", "RUNX1.RUNX1T1", "NUP98.NSD1")
for (gene in geneMuts) {
  formula <- as.formula( paste("Surv(time, status) ~ cyto + wbc.over.16 + strata(age.over.60) + ", gene) ) 
  print( coxph(formula, data = model.data) )
}

# Final Model
coxph(Surv(time, status == 1) ~ cyto + wbc.over.16 + TP53 + FLT3 + strata(age.over.60), data = model.data)



# ** Overall Survival **
# The event looked at is death
# KM
osKM <- survfit(Surv(time2, status2) ~ cyto, type="kaplan-meier", data = model.data)
plot(osKM, main="Overall Survival by Cytogenetic Risk Group", xlab="Overall Survival (Months)", 
     ylab="Survival Probability", col = c(2, 4, 3), lty = c(5, 1, 4), lwd = 2)
legend("topright", legend = c("Intermed", "Good", "Poor"), lty = c(5, 1, 4), col = c(2, 4, 3), lwd = 2)

# Basic Model
coxph(Surv(time2, status2) ~ cyto + wbc.over.16 + strata(age.over.60), data = model.data)

# Add gene terms separately
for (gene in geneMuts) {
  formula <- as.formula( paste("Surv(time2, status2) ~ cyto + wbc.over.16 + strata(age.over.60) + ", gene) ) 
  print( coxph(formula, data = model.data) )
}

# Final Model
coxph(Surv(time2, status2) ~ cyto + wbc.over.16 + TP53 + strata(age.over.60), data = model.data)



# ** Relapse Free Survival **
# Basic Model
coxph(Surv(time, status3 == 1) ~ cyto + wbc.over.16 + strata(age.over.60), data = model.data)

# Add gene terms separately
for (gene in geneMuts) {
  formula <- as.formula( paste("Surv(time, status3 == 1) ~ cyto + wbc.over.16 + strata(age.over.60) + ", gene) ) 
  print( coxph(formula, data = model.data) )
}

# Final Model
coxph(Surv(time, status3 == 1) ~ cyto + wbc.over.16 + strata(age.over.60), data = model.data)



# ** Death-Free (w/o relapse) Survival **
# Basic Model
coxph(Surv(time, status3 == 2) ~ cyto + wbc.over.16 + strata(age.over.60), data = model.data)

# Add gene terms separately
for (gene in geneMuts) {
  formula <- as.formula( paste("Surv(time, status3 == 2) ~ cyto + wbc.over.16 + strata(age.over.60) + ", gene) ) 
  print( coxph(formula, data = model.data) )
}

# Final Model
coxph(Surv(time, status3 == 2) ~ cyto + wbc.over.16 + strata(age.over.60), data = model.data)




