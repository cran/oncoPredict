## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------

library(oncoPredict)

#This script provides an example of how to use calcPhenotype() for drug response prediction as well as its optional #parameters. 

#Set the seed for reproducibility. 
set.seed(12345)

#Determine parameters for calcPhenotype() function.

#Read training data for GDSC (expression and response)
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
trainingPtype = readRDS(file = "GDSC2_Res.rds")
#dim(trainingPtype) #805 198

#GDSC2 expression data for the vignette (it's a much smaller sampling)
trainingExprData=readRDS(file='GDSC2_Expr_short.rds')
#dim(trainingExprData) #1000 400

#IMPORTANT note: here I do e^IC50 since the IC50s are actual ln values/log transformed already, and the calcPhenotype function Paul #has will do a power transformation (I assumed it would be better to not have both transformations)
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
testExprData=as.matrix(read.table('prostate_test_data.txt', header=TRUE, row.names=1))
#dim(testExprData) #20530 550

#Additional parameters. 
#_______________________________________________________
#batchCorrect options: "eb" for ComBat, "qn" for quantiles normalization, "standardize", or "none"
#"eb" is good to use when you use microarray training data to build models on microarray testing data.
#"standardize is good to use when you use microarray training data to build models on RNA-seq testing data (this is what Paul used in the 2017 IDWAS paper that used GDSC microarray to impute in TCGA RNA-Seq data, see methods section of that paper for rationale)
batchCorrect<-"eb"

#Determine whether or not to power transform the phenotype data.
#Default is TRUE.
powerTransformPhenotype<-TRUE

#Determine percentage of low varying genes to remove.
#Default is 0.2 (seemingly arbitrary).
removeLowVaryingGenes<-0.2

#Determine method to remove low varying genes.
#Options are 'homogenizeData' and 'rawData'
#homogenizeData is likely better if there is ComBat batch correction, raw data was used in the 2017 IDWAS paper that used GDSC microarray to impute in TCGA RNA-Seq data.
removeLowVaringGenesFrom<-"homogenizeData"

#Determine the minimum number of training samples required to train on.
#Note: this shouldn't be an issue if you train using GDSC or CTRP because there are many samples in both training datasets.
#10, I believe, is arbitrary and testing could be done to get a better number.
minNumSamples=10

#Determine how you would like to deal with duplicate gene IDs.
#Sometimes based on how you clean the data, there shouldn't be any duplicates to deal with.
#Options are -1 for ask user, 1 for summarize by mean, and 2 for disregard duplicates
selection<- 1

#Determine if you'd like to print outputs.
#Default is TRUE.
printOutput=TRUE

#Indicate whether or not you'd like to use PCA for feature/gene reduction. Options are 'TRUE' and 'FALSE'.
#Note: If you indicate 'report_pca=TRUE' you need to also indicate 'pca=TRUE'
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

#Indicate percent variability (of the training data) you'd like principal components to reflect if pcr=TRUE. Default is .80
percent=80

#Run the calcPhenotype() function using the parameters you specified above.
#__________________________________________________________________________________________________________________________________
wd<-tempdir()
savedir<-setwd(wd)

calcPhenotype(trainingExprData=trainingExprData,
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

#If pcr is performed, you can view a drug's first two principal components (and so on) using the code below. 
#View(load('./calcPhenotype_Output/Vinblastine_1004.RData'))
#View(pcs[,1,1]) #The first pc. 
#View(pcs[,1,2]) #The second pc. 

