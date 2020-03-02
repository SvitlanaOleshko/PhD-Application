rm(list=ls())

setwd("C:/Users/SVITLANA/Desktop/PhD Application/Assignment")

#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("beadarray")
#install.packages("plyr")


#library(beadarray)
load("expression_data.rda")
pheno_data <- read.table("pheno_dat.txt")

#shaping the data
colnames(pheno_data) <- sapply(pheno_data[1,], as.character)
pheno_data <- pheno_data[-1,]

#understand the expression data

#BSData$sampleID 
# sampleID in expression is the same as ID in pheno

head(fData(BSData))
fData(BSData) # features array probes 
table(fData(BSData)[,"Status"]) # only regular

pData(BSData) 
ncol(pData(BSData)) # = 10
length(pData(BSData)[,1]) # 120 => samples

colnames(pData(BSData))

##channelNames(BSData) # only "G"



##################################################################################

# 1.a
assay.Data <- attr(BSData,"assayData")
filt <- apply(assay.Data$Detection < 0.05, 1, sum) >=  dim(assay.Data$Detection)[2]*0.5
BSData.filt <- BSData[filt,]
dim(BSData)
dim(BSData.filt)


#1.b
BSData.norm <- normaliseIllumina(BSData.filt, method="quantile", transform="log2")
dim(BSData.norm)

#1.c
temp <- merge(pData(BSData.norm), 
              pheno_data, by.x = "sampleID", by.y = "ID")
pData(BSData.norm) <- temp

missing <- apply(is.na(pData(BSData.norm)),1,sum) > 0
pData(BSData.norm)[which(missing),"individual_ID"]

BSData.norm.clean <- BSData.norm[,!missing]    

assay.Data.clean <- attr(BSData.norm.clean,"assayData")
filt.clean <- apply(assay.Data.clean$Detection < 0.05, 1, sum) >=  dim(assay.Data.clean$Detection)[2]*0.5
BSData.filt.clean <- BSData.norm.clean[filt.clean,]
dim(BSData.norm.clean)
dim(BSData.filt.clean)

##### ********** CORRECTION *****************************
#library(limma)

#assay.Data.filt.clean$exprs[1:5,1:5]
#assay.Data.filt.clean.batch1Corrected[1:5,1:5]
assay.Data.filt.clean <- attr(BSData.filt.clean,"assayData")
assay.Data.filt.clean.batch1Corrected <- removeBatchEffect(assay.Data.filt.clean$exprs, 
                                          batch = BSData.filt.clean$batch1)
assay.Data.filt.clean.batch12Corrected <- 
  removeBatchEffect(assay.Data.filt.clean.batch1Corrected, 
                    batch = BSData.filt.clean$batch2)

assay.Data.filt.clean.batch123Corrected <- 
  removeBatchEffect(assay.Data.filt.clean.batch12Corrected, 
                    batch = BSData.filt.clean$batch3)

assay.Data.corrected <- assay.Data.filt.clean.batch123Corrected 
assay.Data.corrected[1:5,1:5]

#quality checks

#dim(assay.Data.filt.clean$exprs)
#dim(assay.Data.filt.clean.batch123Corrected)

#assay.Data.filt.clean$exprs[1:5,1:5]
#assay.Data.filt.clean.batch123Corrected[1:5,1:5]

#assay.Data.filt.clean.corrected <- attr(BSData.filt.clean,"assayData")
#assay.Data.filt.clean.corrected$exprs[1:5,1:5]



#BSData.filt.clean.corrected <- 
#              ExpressionSet(assayData, 
#              phenoData=annotatedDataFrameFrom(assayData, byrow=FALSE), 
#              featureData=annotatedDataFrameFrom(assayData, byrow=TRUE), 
#              experimentData=MIAME(), 
#              annotation=character(), 
#              protocolData=annotatedDataFrameFrom(assayData, byrow=FALSE), ...)

#ExpressionSet(assayData = assay.Data.filt.clean.batch123Corrected,
#                      phenoData = pData(BSData.filt.clean),
#                      featureData = fData(BSData.filt.clean))

#BSData.filt.clean.corrected <- assayDataElementReplace(BSData.filt.clean, 
#                                                       "exprs", 
#                                                       assay.Data.filt.clean.batch123Corrected)
#assayData(BSData.filt.clean)

#exampleSummaryData.med = assayDataElementReplace(exampleSummaryData, 
#                                                 "exprs", 
#                                                 exampleSummaryData)
#data(exampleSummaryData)
################################################################################

# 2 - outlier filtering for covariates

## try leverage

treatment <- BSData.filt.clean$treatment
sex <- BSData.filt.clean$sex
age <- BSData.filt.clean$age
bmi <- BSData.filt.clean$bmi


age_num <- as.numeric(as.character(age))
sex_num <- as.numeric(as.character(sex))
treatment_num <- as.numeric(as.character(treatment))
bmi_num <- as.numeric(as.character(bmi))

hist(age_num, 
     col = "lightgreen",
     main = "Histogram of Age Variable",
     xlab = "Age")

hist(bmi_num, 
     col = "lightgreen",
     main = "Histogram of Age Variable",
     xlab = "Age")


X <- as.matrix(data.frame(age_num, bmi_num))

Xt <- t(X)
W <- Xt %*% X
det(W)!=0

# library (MASS)

W_inv <- ginv(W)
H <- X %*% W_inv %*% Xt
leverage <- ifelse(diag(H)>4/118,"high","not high")
# 4/nrow(covariates)


data.temp <- data.frame(age_num, bmi_num,exps1 = as.numeric(t(assay.Data.corrected[1,])))

ggplot(data = data.temp, aes(x=age_num, y=exps1,color=leverage)) +
  geom_point() +
  geom_smooth(data=data.temp,aes(x=age_num, y=exps1),
   method = "lm",inherit.aes = FALSE)

# => leverage captures only covariates outliers (without accounting the association 
# with some outcome variable). However, this sould be enough for us. 

## individual 4 has bmi=200 and age=200 as well as leverage is high

outlier <- leverage == "high"

BSData.final <- BSData.filt.clean[,!outlier]
pData(BSData.final)

## remove outliers from assayData:
assay.Data.corrected <- assay.Data.corrected[,!outlier]
dim(assay.Data.corrected)


### ANOVA
#assay.Data.corrected
#class(assay.Data.corrected)
intensity.trans <- t(assay.Data.corrected)

n.probes <- dim(intensity.trans)[2]
n.samples <- dim(intensity.trans)[1]

temp.mat <- matrix(0, nrow = n.samples, ncol = n.probes)
for (k in 1:n.probes){
  linear.model.fit <- lm(formula = intensity.trans[,k] ~ 
                           treatment_num + sex_num + age_num + bmi_num)
  linear.model.resid <- linear.model.fit$residuals
  resid.mean <- mean(linear.model.resid)
  resid.sd <- sd(linear.model.resid)
  temp.mat[,k] <- abs(linear.model.resid - resid.mean)/resid.sd > 3
  
}

outlier <- apply(temp.mat, 1, sum)/n.probes > 0.05

BSData.final <- BSData.filt.clean[,!outlier]
pData(BSData.final)

## the indiv with age 200 was not removed !!!

###################################################################################

## 3

BSData.rna <- factor(pData(BSData.final)[,"treatment"])
BSData.design <- model.matrix(~0+BSData.rna)
colnames(BSData.design) <- c("Before", "After")
#colnames(BSData.design) <- levels(BSData.rna)
BSData.fit <- lmFit(assay.Data.corrected, BSData.design)
BSData.contrasts <- makeContrasts(Before-After, levels=BSData.design)
BSData.contr.fit <- eBayes(contrasts.fit(BSData.fit, BSData.contrasts))
topTable(BSData.contr.fit, coef=1)



D <- data.frame(fData(BSData.final), pval = BSData.contr.fit$p.value )
D <- D[,c("ProbeID", "SYMBOL", "Before...After")]
names(D) <- c("ProbeID", "SYMBOL", "pval")
sum(D$pval < 0.05) # 138
signif.Genes <- D[D$pval < 0.05,]
unique(signif.Genes$SYMBOL) # 132

######################################################################
#4

# Next, we estimate the fold-change of a gene measured
# by a probe based on the quantile-normalized data. 
# The weighting factor for a probe is calculated based on a Gaussian window function

topTable(BSData.contr.fit, coef=1)
genes <- fData(BSData.final)[,"SYMBOL"]
topTable(BSData.contr.fit, coef=1, number=10000, genelist=genes, 
         adjust.method="none", p.value=0.05, lfc=1.2)

