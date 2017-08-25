# Categorical Method

**Installation**
```
	R CMD SHLIB McCOIL_categorical_code.c llfunction_het.c
```
**Input format**
>Assuming there are n individuals and k loci, SNP calling information is stored in a matrix,
![SNP Calls Matrix](https://cdn.rawgit.com/Greenhouse-Lab/THEREALMcCOIL/gh-pages/assets/equation.svg)
>
>where each element Sij represents SNP information at locus j of individual i, and can be 0 [homozygous minor allele], 0.5 [heterozygous], 1 [homozygous major allele] or -1 [missing data].


**Run**


`See test_R_code.R for an example.`


**Usage**


`McCOIL_categorical(data, ...)`


**Arguments**

|Argument|Description|
|--------|-----------|
data|An R data frame of SNP calling information. Row names are names of samples and column names are names of assays.
maxCOI|Upper bound for COI. The default is 25.
threshold_ind|The minimum number of sites for a sample to be considered. The default is 20.
threshold_site|The minimum number of samples for a locus to be considered. The default is 20.
totalrun|The total number of MCMC iterations. The default is 10000.
burnin|The total number of burnin iterations. The default is 1000.
M0|Initial COI. The default is 15.
e1|The probability of calling homozygous loci heterozygous. The default is 0.05.
e2|The probability of calling heterozygous loci homozygous. The default is 0.05.
path|The default is the current directory.
output|The name of output file. The default is “output.txt”.
err_method|The default is 1.<br>1: use pre-specified e1 and e2 and treat them as constants<br>3: e1 and e2 are estimated with COI and allele frequencies
