## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(oncoPredict)

#This script provides an example of how to control for GLDS in pre-clinical biomarker discovery. 
#Specifically, this script applies glds functions to GDSCv2 data to obtain p-values and beta values for each #drug-gene association. 
#Controlling for GLDS is important because variability in GLDS is evient in cancer cell lines, and controlling for #this variability improves cancer biomarker discovery.

#Set parameters of completeMatrix()
#_____________________________________________________________________
#nPerms=50

#trainingPtype = readRDS(file = "GDSC2_Res.rds")
#There are some NA values so prcomp() will complain when you apply the GLDS function.
#senMat=trainingPtype

#Apply completeMatrix()
#_____________________________________________________________________
#This function will create a file called complete_matrix_output.txt in your working directory.
#This file is used as input in the next function.

#completeMatrix(trainingPtype)

#Apply the glds function.
#_____________________________________________________________________

#Set parameters...

#'@param drugMat A matrix of drug sensitivity data. rownames() are pre-clinical samples, and colnames() are drug names.
#The sensitivity data used here is GDSCv2.

#Read GDSC's updated cell line information file (used later).
#cellLineDetails<-read_excel('Cell_Lines_Details.xlsx')
cellLineDetails<-read.csv('Cell_Lines_Details.csv')

#I ran the response data through completeMatrix() because there were some NA values.
#NA values in the response data will cause a problem when you apply prcomp().
cm<-read.table('complete_matrix_output_GDSCv2.txt', header=TRUE, row.names=1) #Now, there are no NA values.

#rownames(cm) #Cosmic identifiers are used for cell names in this dataset...this will cause a problem later when matching cell lines between sensitivity and mutation data.

#Replace the rownames of cm with cell line names. Right now, they are cosmic ids.
#This will require using GDSC's cell line details file (which maps cosmic ids to cell line names).
newRows <- substring(rownames(cm),8) #Remove 'COSMIC'...keep the numbers after COSMIC.
indices<-match(as.numeric(newRows), as.vector(unlist(cellLineDetails[,2]))) #Refer to the cell line details file to make this replacement.
newNames<-as.vector(unlist(cellLineDetails[,1]))[indices] #Reports the corresponding cell line names
rownames(cm)<-newNames

#Fix the drug names in the cm object so that it's just the name of the drug (remove those extra numbers/identifiers at the end).
#gdscv2_drugs.xlsx contains the colnames of cm in the correct order, but with the extra identifiers removed.
#fix<-read_excel('gdscv2_drugs.xlsx')
#fix<-as.vector(unlist(fix[,2]))
fix<-as.vector(unlist(read.table('gdscv2_drugs.txt', header=TRUE)))
colnames(cm)<-as.vector(fix)
drugMat<-as.matrix(cm) #Finally, set this object as the drugMat parameter. 
#dim(drugMat) #805 samples vs. 198 drugs

#'@param markerMat A matrix containing the data for which you are looking for an association with drug sensitivity (e.g. a matrix of somatic mutation data). rownames() are 
#marker names (e.g. gene names), and colnames() are samples.
#The dataset used here is GDSCv2's updated mutation data for pan-cancer. It includes both CNV and coding variant.
#mutationMat<-read.csv('GDSC2_Pan_Both.csv')
#mutationMat<-mutationMat[,c(1,6,7)] #Index to these 3 columns of interest.
#colnames(mutationMat) #"cell_line_name"  "genetic_feature" "is_mutated" 
#Some entries are duplicated cell line name - genetic feature combos...remove them to avoid problems with pivot_wider().
#vec<-c()
#for (i in 1:nrow(mutationMat)){
#  vec[i]<-paste(mutationMat[i,1],mutationMat[i,2], sep=' ')
#}
#nonDupIndices<-match(unique(vec), vec)
#mutationMat2<-mutationMat[nonDupIndices,]

#Some gene mutation entries are blank...remove them to avoid problems with pivot_wider().
#library(tidyverse)
#good<-(mutationMat2[,2]) != ""
#mutationMat3<-mutationMat2[good,]
#mutationMat4<-mutationMat3 %>%
#  pivot_wider(names_from=genetic_feature,
#              values_from=is_mutated)
#rownames(mutationMat4)<-as.vector(unlist(mutationMat4[,1])) #Make cell lines the #rownames...right now they are column 1.
#cols<-rownames(mutationMat4)
#mutationMat4<-as.matrix(t(mutationMat4[,-1]))
#Make sure the matrix is numeric.
#mutationMat<-mutationMat4
#mutationMat4<-apply(mutationMat4, 2, as.numeric)
#rownames(mutationMat4)<-rownames(mutationMat)
#markerMat<-mutationMat4
# replace all non-finite values with 0
#markerMat[!is.finite(markerMat)] <- 0
#colnames(markerMat)<-cols
#dim(markerMat) #1315 1389
#write.table(markerMat, file='markerMat.txt')
markerMat<-as.matrix(read.table('markerMat.txt', header=TRUE, row.names=1))

#'@param drugRelatedness 
#This file is GDSC's updated drug relatedness file (obtained from bulk data download/all compounds screened/compounds-annotation).
#Note: I had to change some drug names in this file so that they matched colnames of cm.
#Ex: replace - with . (small modifications like that).
drugRelatedness <- read.csv("screened_compunds_rel_8.2.csv")
drugRelatedness<-drugRelatedness[,c(3,6)]
#colnames(drugRelatedness) #"DRUG_NAME"      "TARGET_PATHWAY"

wd<-tempdir()
savedir<-setwd(wd)

glds(drugMat,
     drugRelatedness,
     markerMat,
     minMuts=5,
     additionalCovariateMatrix=NULL,
     threshold=0.7)

