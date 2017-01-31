path= "../proportional_method"  ##enter your path here
setwd(path)
source("McCOIL_proportional.R")

#example dataset
dataA1i = read.table(paste(path, "/dataA1_test.txt", sep=""), head=T)
dataA2i = read.table(paste(path, "/dataA2_test.txt", sep=""), head=T)
dataA1= dataA1i[,-1]
dataA2= dataA2i[,-1]
rownames(dataA1)= dataA1i[,1]
rownames(dataA2)= dataA2i[,1]

McCOIL_proportional(dataA1, dataA2, maxCOI=25, totalrun=5000, burnin=100, M0=15, epsilon=0.02, err_method=3, path=getwd(), output="output_test.txt" )

