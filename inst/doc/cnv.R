## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------

library(oncoPredict)

#This script provides an example of how to download cnv (copy number variation) data from the GDC database for GBM #(glioblastoma), how to apply the map_cnv() function to that data to map cnvs to genes, and how to test the drugs in your drug #response dataset to each cnv to identify biomarkers that enrich for drug response.

#First, download CNV data for your cancer of interest from GDC database. The cnv data will be exported to your #working directory as cnv.txt 

#This code will export the CNV data into a text file called, 'cnv.txt', containing a table with colnames() 'Sample', 'Chromosome', 'Start', 'End', 'Num Probes', 'Segment_Mean'
#The genome of reference is hg19.
#query.gbm.nocnv<-GDCquery(project = "TCGA-GBM",
#                          data.category = "Copy number variation",
#                          legacy = TRUE,
#                          file.type = "nocnv_hg19.seg",
#                          sample.type = c("Primary Tumor"))
#patient_total<-nrow((query.gbm.nocnv$results)[[1]]) #The total number of patients GDC has CNV data for
#query.gbm.nocnv$results[[1]]<-query.gbm.nocnv$results[[1]][1:patient_total,]
#GDCdownload(query.gbm.nocnv, files.per.chunk = 100)
#gbm.nocnv<-GDCprepare(query.gbm.nocnv, save = TRUE, save.filename = "GBMnocnvhg19.rda")
#write.table(gbm.nocnv, file='cnv.txt')

#Second, apply map_cnv() to map cnv data to genes. The mapped cnv data will be exported to your working directory as map.RData
#The mapping is accomplished by intersecting the gene with the overlapping CNV level. If the gene isn't fully #captured by the CNV, an NA will be assigned.

#Determine the parameters of the map_cnv() function.
Cnvs<-read.table('cnv.txt', header=TRUE, row.names=1)

#Third, apply idwas() to test each cnv and each drug. The p-values and beta-values for each test will be exported to #your working directory as CnvTestOutput_pVals.txt and CnvTestOutput_betas.txt

#Determine the parameters of the idwas() function...

#Set the drug_prediction parameter.
#Make sure rownames() are samples, and colnames() are drugs. Also make sure this data is a data frame.
drug_prediction<-t(as.data.frame(read.table('DrugPredictions.txt', header=TRUE, row.names=1)))
#dim(drug_prediction) #165 198

#In this example, I had to replace the '.' in the names of these TCGA samples with '-' so that they are of the same form as samples in the cnv  data (you may not have to do this).
rownames(drug_prediction)<-gsub(".", "-", rownames(drug_prediction), fixed=T)

#Make sure the sample identifiers in the 'drug prediction' data are of similar form as the sample identifiers in the 'data' parameter.
rows=rownames(drug_prediction)
rownames(drug_prediction)<-substring(rows, 3, nchar(rows))
drug_prediction<-as.data.frame(drug_prediction)

#Determine the number of samples you want the CNVs to be amplified in. The default is 10.
n=10

#Indicate whether or not you would like to test cnv data. If TRUE, you will test cnv data. If FALSE, you will test mutation data.
cnv=TRUE

wd<-tempdir()
savedir<-setwd(wd)

#Apply map_cnv()
#This function produces the file map.RData, which stores the object 'theCnvQuantVecList_mat'
map_cnv(Cnvs=Cnvs)

#Set the data parameter.
load('map.RData') #This loads the object 'theCnvQuantVecList_mat', which was obtained using map_cnv()
#Make sure this data is a data frame and that colnames() are samples.
data<-as.data.frame(t(theCnvQuantVecList_mat))
samps<-colnames(data)
colnames(data)<-substr(samps,1,nchar(samps)-12)

#Apply idwas()
idwas(drug_prediction=drug_prediction, data=data, n=n, cnv=cnv)


