rm(list=ls())

setwd("C:/Users/SVITLANA/Desktop/PhD Application/Assignment")

#if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
#BiocManager::install("beadarray")

library(beadarray)
library(limma)
library(MASS)

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
channelNames(BSData) # only "G"

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

###### ***************************** CORRECTION ***************************** ######

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

################################################################################

# 2 - outlier filtering for covariates

## leverage

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

# => leverage captures only covariates outliers 
# (without accounting the association with some outcome variable). 
# However, this should be enough for us. 

# individual 4 has bmi=200 and age=200 as well as leverage is high

outlier <- leverage == "high"

BSData.final <- BSData.filt.clean[,!outlier]
pData(BSData.final)

## remove outliers from assayData:
assay.Data.corrected <- assay.Data.corrected[,!outlier]
dim(assay.Data.corrected)


## ANOVA

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

## the indiv with age 200 was not removed => stay with leverage method

###################################################################################

## 3

rna <- factor(pData(BSData.final)[,"treatment"])
# the same as BSData.final$treatment
# However, rna has levels "0", "1"
# when BSData.final$treatment has levels "0", "1", "treatment"

design <- model.matrix(~0+rna)
# matrix 112x2: 1st column is "treatment = F", 2nd column is "treatment = T"

colnames(design) <- c("Before", "After")

aw <- arrayWeights(assay.Data.corrected, design)
# Estimates relative quality weights for each array in a multi-array experiment.

fit <- lmFit(assay.Data.corrected, design, weights=aw)
contrasts <- makeContrasts(Before-After, levels=design)
contr.fit <- eBayes(contrasts.fit(fit, contrasts))
topTable(contr.fit, coef=1)


D <- data.frame(fData(BSData.final), pval = contr.fit$p.value )
D <- D[,c("ProbeID", "SYMBOL", "Before...After")]
names(D) <- c("ProbeID", "SYMBOL", "pval")
sum(D$pval < 0.05) 
# 335
signif.Genes <- D[D$pval < 0.05,]
unique(signif.Genes$SYMBOL) 
# 307

####################################################################################
#4

# Next, we estimate the fold-change of a gene measured
# by a probe based on the quantile-normalized data. 
# The weighting factor for a probe is calculated based on a Gaussian window function

topTable(contr.fit, coef=1)
genes <- fData(BSData.final)[,"SYMBOL"]
topTable(contr.fit, coef=1, number=10000, genelist=genes, 
         adjust.method="none", p.value=0.05, lfc=log2(1.2))

log2(1.2) # =0.2630344

