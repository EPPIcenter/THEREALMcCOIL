McCOIL_proportional = function(dataA1, dataA2, maxCOI=25, totalrun=10000, burnin=1000, M0=15, epsilon=0.02, path=getwd(), output="output.txt" ){

	grid = read.table(paste(path, "/fitted_beta_grid_25.txt", sep=""), head=T)
	n=nrow(dataA1)
	k=ncol(dataA1)
	M0=rep(M0, n)
	P0=rep(0.5, k)
	A1=as.vector(t(dataA1))
	A2=as.vector(t(dataA2))

	if ((n>10 & k>10)){	
		dyn.load(paste(path, "/McCOIL_prop_code.so", sep=""))
		Kc <- .C("McCOIL_prop", as.integer(maxCOI), as.integer(totalrun), as.integer(n), as.integer(k), as.double(A1), as.double(A2), as.integer(M0), as.double(P0), as.double(grid$A), as.double(grid$B), as.double(epsilon), as.character(output), as.character(path))
		dyn.unload(paste(path, "/McCOIL_prop_code.so", sep=""))

	} else { stop(paste("Sample size is too small (n=", n, ", k=", k,").", sep=""))}
		
	##summarize results
	outputMCMC2 = read.table(paste(path, "/", output, sep=""), head=F)
	meanM= as.numeric(round(apply(outputMCMC2[(burnin+1): totalrun, (1:n)+1], 2, mean)))
	meanP= as.numeric(apply(outputMCMC2[(burnin+1): totalrun, ((1:k)+n+1)], 2, mean))
	medianM= as.numeric(apply(outputMCMC2[(burnin+1): totalrun, (1:n)+1], 2, median))
	medianP= as.numeric(apply(outputMCMC2[(burnin+1): totalrun, ((1:k)+n+1)], 2, median))
	M975= as.numeric(apply(outputMCMC2[(burnin+1): totalrun, (1:n)+1], 2, function(x) quantile(x, probs= 0.975)))
	P975= as.numeric(apply(outputMCMC2[(burnin+1): totalrun, ((1:k)+n+1)], 2, function(x) quantile(x, probs= 0.975)))
	M025= as.numeric(apply(outputMCMC2[(burnin+1): totalrun, (1:n)+1], 2, function(x) quantile(x, probs= 0.025)))
	P025= as.numeric(apply(outputMCMC2[(burnin+1): totalrun, ((1:k)+n+1)], 2, function(x) quantile(x, probs= 0.025)))
	sdM= as.numeric(apply(outputMCMC2[(burnin+1): totalrun, (1:n)+1], 2, sd))
	sdP= as.numeric(apply(outputMCMC2[(burnin+1): totalrun, ((1:k)+n+1)], 2, sd))

	output_sum= data.frame(cbind(rep(output, (n+k)), c(rep("C", n), rep("P", k)), c(rownames(dataA1), colnames(dataA1)), c(meanM, meanP), c(medianM, medianP), round(c(sdM, sdP), digits=5), c(M025, P025), c(M975, P975)))
	colnames(output_sum)=  c("file", "CorP","name","mean","median","sd", "quantile0.25", "quantile0.975")
	write.table(output_sum, paste(path, "/", output, "_summary.txt", sep=""), sep="\t", col.names=T, row.names=F, quote=F)

}


