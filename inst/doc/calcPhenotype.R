## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

vignette_file <- function(...) {
  candidates <- c(
    file.path(...),
    file.path("vignettes", ...),
    file.path("inst", "extdata", ...),
    file.path(Sys.getenv("PWD"), "inst", "extdata", ...),
    system.file("extdata", ..., package = "oncoPredict"),
    system.file("doc", ..., package = "oncoPredict")
  )
  candidates <- candidates[nzchar(candidates) & file.exists(candidates)]
  if (!length(candidates)) {
    stop("Could not find vignette file: ", file.path(...), call. = FALSE)
  }
  candidates[[1]]
}

## ----setup--------------------------------------------------------------------

library(oncoPredict)

#This vignette demonstrates how to use calcPhenotype() for drug response prediction and how to set commonly used optional parameters.

#Set the seed for reproducibility.
set.seed(12345)

#Determine parameters for calcPhenotype() function.

#Read training data for GDSC (expression and response).
#_______________________________________________________

#GDSC1 Data
#Read GDSC training expression data. rownames() are genes and colnames() are samples (cell lines/cosmic ids).
#trainingExprData=(readRDS("GDSC1_Expr.rds"))
#dim(trainingExprData) #17419 958
#Read GDSC1 response data. rownames() are samples (cell lines, cosmic ids), colnames() are drugs.
#trainingPtype = (readRDS("GDSC1_Res.rds"))
#dim(trainingPtype) #958 367 For GDSC1
#trainingPtype<-trainingPtype[,1:2] #Just 2 drugs for the vignette.

#GDSC2 Data
#Read GDSC training expression data. rownames() are genes and colnames() are samples.
#trainingExprData=readRDS(file='GDSC2_Expr.rds')
#dim(trainingExprData) #17419 805
#Read GDSC2 response data. rownames() are samples, colnames() are drugs.
trainingPtype = readRDS(file = vignette_file("GDSC2_Res.rds"))
#dim(trainingPtype) #805 198

#GDSC2 expression data for this vignette. This is a reduced example dataset.
trainingExprData=readRDS(file = vignette_file("GDSC2_Expr_short.rds"))
#dim(trainingExprData) #500 200

#The GDSC IC50 values are already log-transformed. Convert them back before using
#the default power transformation in calcPhenotype().
trainingPtype<-exp(trainingPtype)

#Or read training data for CTRP (expression and response)
#_______________________________________________________
#Read CTRP training expression data. rownames() are genes and colnames() are samples (cell lines/cosmic ids).
#trainingExprData = readRDS(file = "CTRP2_Expr.rds")
#dim(trainingExprData) #51847 829
#Read CTRP training response data. rownames() are samples (cell lines, cosmic ids), colnames() are drugs.
#trainingPtype = readRDS(file = "CTRP2_Res.rds")
#dim(trainingPtype) #829 545

#Test data.
#_______________________________________________________
#Read testing data as a matrix with rownames() as genes and colnames() as samples.
testExprData=as.matrix(read.table(vignette_file("prostate_test_data.txt"), header=TRUE, row.names=1))
#dim(testExprData) #1000 20

#Additional parameters.
#_______________________________________________________
#batchCorrect options: "eb" for ComBat, "qn" for quantile normalization, "standardize" for z-score standardization, "rank", "rank_then_eb", or "none"
#"eb" is often used when the training and testing datasets are both microarray data.
#"standardize" can be useful when training with microarray data and predicting in RNA-seq data.
batchCorrect<-"eb"

#Determine whether or not to power transform the phenotype data.
#Default is TRUE.
powerTransformPhenotype<-TRUE

#Determine percentage of low varying genes to remove.
#Default is 0.2. This filter reduces the influence of low-variance genes and can
#also reduce runtime for large RNA-seq matrices.
removeLowVaryingGenes<-0.2

#Determine method to remove low varying genes.
#Options are 'homogenizeData' and 'rawData'
#Use 'homogenizeData' to filter after the training and test matrices have been
#aligned and batch-corrected.
removeLowVaringGenesFrom<-"homogenizeData"

#Determine the minimum number of training samples required to train on.
#This example uses reduced data, so the threshold is set low enough for the
#vignette. Larger analyses should use enough samples to fit a reliable model.
minNumSamples=10

#Determine how you would like to deal with duplicate gene IDs.
#Depending on preprocessing, duplicate gene identifiers may or may not be present.
#Options are -1 for ask user, 1 for summarize by mean, and 2 for disregard duplicates
selection<- 1

#Determine if you'd like to print outputs. Set to FALSE here to keep the vignette output concise.
printOutput=FALSE

#Indicate whether or not you'd like to use principal component regression for feature/gene reduction. Options are 'TRUE' and 'FALSE'.
#Note: If you indicate 'report_pc=TRUE' you need to also indicate 'pcr=TRUE'
pcr=FALSE

#Indicate whether you want to output the principal components. Options are 'TRUE' and 'FALSE'.
report_pc=FALSE

#Indicate if you want correlation coefficients for biomarker discovery. These are the correlations between a given gene of interest across all samples vs. a given drug response across samples.
#These correlations can be ranked to obtain a ranked correlation to determine highly correlated drug-gene associations.
cc=FALSE

#Indicate whether or not you want to output the R^2 values for the data you train on from true and predicted values.
#These values represent the percentage in which the optimal model accounts for the variance in the training data.
#Options are 'TRUE' and 'FALSE'.
rsq=FALSE

#Indicate percent variability (of the training data) you'd like principal components to reflect if pcr=TRUE. Default is 80.
percent=80

#Run calcPhenotype() using the parameters specified above.
#__________________________________________________________________________________________________________________________________
drug_predictions <- calcPhenotype(trainingExprData=trainingExprData,
                                  trainingPtype=trainingPtype,
                                  testExprData=testExprData,
                                  batchCorrect=batchCorrect,
                                  powerTransformPhenotype=powerTransformPhenotype,
                                  removeLowVaryingGenes=removeLowVaryingGenes,
                                  minNumSamples=minNumSamples,
                                  selection=selection,
                                  printOutput=printOutput,
                                  pcr=pcr,
                                  removeLowVaringGenesFrom=removeLowVaringGenesFrom,
                                  report_pc=report_pc,
                                  cc=cc,
                                  percent=percent,
                                  rsq=rsq)

#The returned matrix contains predicted drug response values for each test sample.
#dim(drug_predictions)

#Visualize the predicted drug response matrix.
drug_predictions_scaled <- t(scale(drug_predictions))
drug_predictions_scaled[!is.finite(drug_predictions_scaled)] <- 0
heatmap(drug_predictions_scaled,
        Rowv=NA,
        Colv=NA,
        scale="none",
        labRow=NA,
        cexCol=0.7,
        margins=c(6, 2),
        xlab="Test samples",
        ylab="Drugs",
        main="Predicted drug response")

