#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include "loglikelihood_het.h"
double logLike_het(int M, double P, double S2, double e1, double e2){
	double lltrue = 0.0;
	if (S2 ==1){
		lltrue= log( (1-e1)*((double)pow(P,M)) + e2/2*(1- (double)pow(P,M)- (double)pow((1-P),M)) );
	}
	else if (S2 ==0){
		lltrue= log( (1-e1)*((double)pow((1-P),M)) + e2/2*(1- (double)pow(P,M)- (double)pow((1-P),M)) );
	}
	else if (S2 ==0.5){
		if ((M==1) && (e1==0)) lltrue=-9999999;
		else lltrue= log( e1*((double)pow(P,M)) + e1*((double)pow((1-P),M)) + (1-e2)*(1- (double)pow(P,M)- (double)pow((1-P),M)) );
	}
	else lltrue=0.0;
	return (lltrue);


}

