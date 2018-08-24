# THEREALMcCOIL

THE REAL McCOIL is a Markov chain Monte Carlo method that estimates complexity of infection and population allele frequencies using SNP data obtained from Sequenom or similar types of SNP assays.  

We used the data in two ways:

1. categorical method, in which we considered SNP to be heterozygous or homozygous

2. proportional method, in which the relative signal intensity for each allele is used.

The codes can be run on both Mac and Windows. After compiling, you can run the code using the function McCOIL_categorical (categorical) or McCOIL_proportional (proportional) in R.
The instructions of how to compile the codes can be found in README.docx in the folders of categorical method and proportional method. You can also find example R codes that show how to use R functions as well as example input and output files in those folders (example R code: test_R_code.R; input files: input_test.txt [categorical], dataA1_test.txt and dataA2_test.txt [proportional]; output files: output_test.txt and output_test.txt_summary.txt). Please note that it is important to check the convergence of MCMC and it is recommended to run the codes multiple times and check for consistency between them.

***
01-31-2017: THE REAL McCOIL v2 includes the function of estimating parameters of measurement error with COI and allele frequencies.
***
03-14-2018: There is a web version of [THEREALMcCOIL](http://35.196.107.63/run). 
***
08-24-2018: Some explanation for the output file.
In the output file "output_test.txt", each row represents one MCMC run, the first column shows which run it is, and the second to the last columns show the values of different parameters for that run in an order as follows:

COI of individual 1, COI of individual 2, ..., COI of last individual, allele frequency of allele 1, allele frequency of allele 2, ...., allele frequency of last allele.

The last row shows the number of accepted changes out of all MCMC runs for each parameter.
