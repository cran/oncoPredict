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

#This vignette demonstrates how to control for general levels of drug sensitivity
#(GLDS) in pre-clinical biomarker discovery. The example applies glds() to GDSC2
#data to obtain p-values and beta values for drug-marker associations.

#Set parameters of completeMatrix().
#_____________________________________________________________________
#nPerms=50

#trainingPtype = readRDS(file = "GDSC2_Res.rds")
#There are some NA values, which will cause prcomp() to fail when applying GLDS.
#senMat=trainingPtype

#Apply completeMatrix()
#_____________________________________________________________________
#This function returns the completed matrix. Set folder=TRUE to also write complete_matrix_output.txt.

#completeMatrix(trainingPtype)

#Apply the glds() function.
#_____________________________________________________________________

#Set parameters...

#drugMat is a matrix of drug sensitivity data. rownames() are pre-clinical samples, and colnames() are drug names.
#The sensitivity data used here is GDSCv2.

#Read GDSC's updated cell line information file (used later).
#cellLineDetails<-read_excel('Cell_Lines_Details.xlsx')
cellLineDetails<-read.csv(vignette_file("Cell_Lines_Details.csv"))

#The response data were processed with completeMatrix() because NA values in the
#response matrix will cause prcomp() to fail.
cm<-read.table(vignette_file("complete_matrix_output_GDSCv2.txt"), header=TRUE, row.names=1) #No NA values remain.

#Cosmic identifiers are used for cell names in this dataset and are converted
#to cell-line names before matching with the marker matrix.

#Replace the rownames of cm with cell line names. Right now, they are cosmic ids.
#This will require using GDSC's cell line details file (which maps cosmic ids to cell line names).
newRows <- substring(rownames(cm),8) #Remove 'COSMIC'...keep the numbers after COSMIC.
indices<-match(as.numeric(newRows), as.vector(unlist(cellLineDetails[,2]))) #Refer to the cell line details file to make this replacement.
newNames<-as.vector(unlist(cellLineDetails[,1]))[indices] #Reports the corresponding cell line names
# Match the sanitized cell-line names used in the example marker matrix.
rownames(cm)<-make.names(newNames)

#Update the drug names in cm by removing extra identifiers appended to the names.
#gdscv2_drugs.xlsx contains the colnames of cm in the correct order with those identifiers removed.
#fix<-read_excel('gdscv2_drugs.xlsx')
#fix<-as.vector(unlist(fix[,2]))
fix<-as.vector(unlist(read.table(vignette_file("gdscv2_drugs.txt"), header=TRUE)))
colnames(cm)<-as.vector(fix)
drugMat<-as.matrix(cm) #Finally, set this object as the drugMat parameter. 
#dim(drugMat) #100 samples vs. 198 drugs in this reduced example file.

#markerMat contains the data to test for association with drug sensitivity (e.g. a matrix of somatic mutation data). rownames() are
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
#rownames(mutationMat4)<-as.vector(unlist(mutationMat4[,1])) #Use cell lines as rownames before transposing.
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
#write.table(markerMat, file='markerMat.txt')
#The included example markerMat is reduced to the top 200 markers among the samples used here.
markerMat<-as.matrix(read.table(vignette_file("markerMat.txt"), header=TRUE, row.names=1, check.names=FALSE))
#dim(markerMat) #200 markers vs. 40 samples in this reduced example file.

if(length(intersect(colnames(markerMat), rownames(drugMat))) == 0){
  stop("No overlapping samples were found between markerMat and drugMat.")
}

#drugRelatedness contains drug names and the corresponding target pathways.
#This file is GDSC's updated drug relatedness file (obtained from bulk data download/all compounds screened/compounds-annotation).
#Some drug names in this file were adjusted so they match colnames of cm.
#Ex: replace - with . (small modifications like that).
drugRelatedness <- read.csv(vignette_file("screened_compunds_rel_8.2.csv"))
drugRelatedness<-drugRelatedness[,c(3,6)]
#colnames(drugRelatedness) #"DRUG_NAME"      "TARGET_PATHWAY"

glds_results <- glds(drugMat,
                     drugRelatedness,
                     markerMat,
                     minMuts=5,
                     additionalCovariateMatrix=NULL,
                     threshold=0.7)

