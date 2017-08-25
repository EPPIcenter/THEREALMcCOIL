# Proportional Method
**Installation**
```
	R CMD SHLIB McCOIL_prop_code.c llfunction.c
```
**Input format**
>Assuming there are n individuals and k loci, SNP information is stored in two n × k matrices A1 and A2, representing the intensity of two alleles in the SNP assay, respectively. Each element in the matrices, A1ij or A2ij , represents the intensity of allele 1 or allele 2 at locus j of individual i, and missing data is set to -1.

**Run**


`See test_R_code.R for an example.`


**Usage**


`McCOIL_proportional(dataA1, dataA2, ...)`


**Arguments**

|Argument|Description|
|--------|-----------|
|dataA1, dataA2| An R data frame. The intensity of signals of allele 1 and allele 2 from the SNP assay. Row names are names of samples and column names are names of assays.|
|maxCOI|Upper bound for COI. The default is 25.|
|totalrun|The total number of MCMC iterations. The default is 10000.|
|burnin|The total number of burnin iterations. The default is 1000.|
|M0|Initial COI. The default is 15.|
|epsilon|The level of measurement error (εest). The default is 0.2.|
|path|The default is the current directory.|
|output|The name of output file. The default is “output.txt”.|
err_method|The default is 1.<br>1: use pre-specified epsilon and treat as constant.<br>3: epsilon is estimated with COI and allele frequencies.|
