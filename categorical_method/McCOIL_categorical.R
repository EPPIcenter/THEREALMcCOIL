McCOIL_categorical = function(data, maxCOI=25, threshold_ind=20, threshold_site=20, totalrun=10000, burnin=1000, M0=15, e1=0.05, e2=0.05, err_method=1, path=getwd(), output="output.txt" ){

	mcCoil_categorical_code_location = '/McCOIL_categorical_code.so'
	if(Sys.info()['sysname'] == 'Windows') {
		mcCoil_categorical_code_location = '/McCOIL_categorical_code.dll'
	}

	In_ind= rep(NA, nrow(data))
	In_site= rep(NA, ncol(data))
	for (i in (1:nrow(data))){
		if ((length(data[i,])-sum(data[i,]==-1)) >= threshold_ind) In_ind[i]= TRUE
		else In_ind[i]= FALSE
	}
	for (i in (1:ncol(data))){
		if ((length(data[,i])-sum(data[,i]==-1)) >= threshold_site) In_site[i]= TRUE
		else In_site[i]= FALSE
	}
	
	##remove sites and individuals with too much missing data
	simpleS = data[In_ind,]
	simpleS = simpleS[,In_site]
	
	##remove sites with P=0 or 1
	P=rep(NA, ncol(simpleS))
	for (j in (1:ncol(simpleS))){
		temp= simpleS[,j]
		P[j]= (sum(temp==1) + 0.5*sum(temp!=0 & temp!=1& temp!=-1))/sum(temp!=-1)
	}
	In = (P!=Inf & P!="NaN" & P!=0 & P!=1)
	simpleS2= simpleS[, In]
	
	select_pos =colnames(data)[In_site][In]
	select_ind = rownames(data)[In_ind]

	n= nrow(simpleS2)
	k= ncol(simpleS2)
	simpleS2_vec= as.vector(t(simpleS2))
	P0= P[In]
	M0=rep(M0, n)

	if ((n>10 & k>10)){	
		dyn.load(paste(path, mcCoil_categorical_code_location, sep=""))
		K <- .C("McCOIL_categorical", as.integer(maxCOI), as.integer(totalrun), as.integer(n), as.integer(k), as.double(simpleS2_vec), as.integer(M0), as.double(P0), as.double(e1), as.double(e2), as.character(output), as.character(path), as.integer(err_method))
		dyn.unload(paste(path, mcCoil_categorical_code_location, sep=""))
	} else { stop(paste("Sample size is too small (n=", n, ", k=", k,").", sep=""))}

	##summarize results
	outputMCMC1 = read.table(paste(path, "/", output, sep=""), head=F)
	meanM= as.numeric(round(apply(outputMCMC1[(burnin+1): totalrun, (1:n)+1], 2, mean)))
	meanP= as.numeric(apply(outputMCMC1[(burnin+1): totalrun, ((1:k)+n+1)], 2, mean))
	medianM= as.numeric(apply(outputMCMC1[(burnin+1): totalrun, (1:n)+1], 2, median))
	medianP= as.numeric(apply(outputMCMC1[(burnin+1): totalrun, ((1:k)+n+1)], 2, median))
	M975= as.numeric(apply(outputMCMC1[(burnin+1): totalrun, (1:n)+1], 2, function(x) quantile(x, probs= 0.975)))
	P975= as.numeric(apply(outputMCMC1[(burnin+1): totalrun, ((1:k)+n+1)], 2, function(x) quantile(x, probs= 0.975)))
	M025= as.numeric(apply(outputMCMC1[(burnin+1): totalrun, (1:n)+1], 2, function(x) quantile(x, probs= 0.025)))
	P025= as.numeric(apply(outputMCMC1[(burnin+1): totalrun, ((1:k)+n+1)], 2, function(x) quantile(x, probs= 0.025)))
	sdM= as.numeric(apply(outputMCMC1[(burnin+1): totalrun, (1:n)+1], 2, sd))
	sdP= as.numeric(apply(outputMCMC1[(burnin+1): totalrun, ((1:k)+n+1)], 2, sd))	

	if (err_method==3){
		mean_e1= as.numeric(mean(outputMCMC1[(burnin+1): totalrun, (k+n+2)]))
		median_e1= as.numeric(median(outputMCMC1[(burnin+1): totalrun, (k+n+2)]))
		e1_975=  as.numeric(quantile(outputMCMC1[(burnin+1): totalrun,  (k+n+2)], probs= 0.975))
		e1_025=  as.numeric(quantile(outputMCMC1[(burnin+1): totalrun,  (k+n+2)], probs= 0.025))
		sd_e1=  as.numeric(sd(outputMCMC1[(burnin+1): totalrun, (k+n+2)]))
	
		mean_e2= as.numeric(mean(outputMCMC1[(burnin+1): totalrun, (k+n+3)]))
		median_e2= as.numeric(median(outputMCMC1[(burnin+1): totalrun, (k+n+3)]))
		e2_975=  as.numeric(quantile(outputMCMC1[(burnin+1): totalrun,  (k+n+3)], probs= 0.975))
		e2_025=  as.numeric(quantile(outputMCMC1[(burnin+1): totalrun,  (k+n+3)], probs= 0.025))
		sd_e2=  as.numeric(sd(outputMCMC1[(burnin+1): totalrun, (k+n+3)]))	
	}
	if ((err_method==1) | (err_method==2)) {
			output_sum= data.frame(cbind(rep(output, (n+k)),	
							c(rep("C", n), rep("P", k)), c(select_ind, select_pos), c(meanM, meanP), c(medianM, medianP), round(c(sdM, sdP),digits=5), c(M025, P025), c(M975, P975)))
	}
	else {
			output_sum= data.frame(cbind(rep(output, (n+k+2)),	
							c(rep("C", n), rep("P", k), "e1", "e2"), c(select_ind, select_pos, "e1","e2"), c(meanM, meanP, mean_e1, mean_e2), 
							c(medianM, medianP, median_e1, median_e2), round(c(sdM, sdP, sd_e1, sd_e2),digits=5), c(M025, P025, e1_025, e2_025), c(M975, P975, e1_975, e2_975)))
	}
	colnames(output_sum)=  c("file", "CorP","name", "mean","median","sd", "quantile0.025", "quantile0.975")
	write.table(output_sum, paste(path, "/", output, "_summary.txt", sep=""), sep="\t", col.names=T, row.names=F, quote=F)

}


