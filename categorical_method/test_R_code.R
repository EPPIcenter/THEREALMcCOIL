path= "../categorical_method" ##enter your path here
setwd(path)
source("McCOIL_categorical.R")

#example dataset
data0= read.table(paste(path, "/input_test.txt", sep=""), head=T)
data=data0[,-1]
rownames(data)=data0[,1]
McCOIL_categorical(data,maxCOI=25, threshold_ind=20, threshold_site=20, totalrun=1000, burnin=100, M0=15, e1=0.05, e2=0.05, path=getwd(), output="output_test.txt" )

