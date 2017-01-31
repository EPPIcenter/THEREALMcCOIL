//R CMD SHLIB McCOIL_prop_code.c llfunction.c

#include <stdio.h>
#include <R.h>
#include <Rmath.h>
#include "loglikelihood.h"
#include <time.h>

void McCOIL_prop (int *max, int *iterations, int *n0, int *k0, double *A1, double *A2, int *M0, double *P0, double *A, double *B, double* c0, char **file_index, char **path, int *err_method0) {

////////can be changed///////
	double varP=0.1;  
	double varS=0.1;
	double varE=0.0001;
	double upper_bound_c=0.1; //epsilon
	
///////don't change////////
	double temp_l_est=0;
	double temp_l_true=0;

	int i=0, j=0, prime=0, x=0, y=0;
	double sumori=0, sumcan=0; 
	
	GetRNGstate();  //provide seed for random sampling
	
	////read in the values////
	double c = *c0;  //epsilon
	int max_moi= *max;
	int iter = *iterations;
	int n = *n0, k = *k0;
	int err_method= *err_method0; //1: use pre-specified epsilon; 2: use likelihood-free sampling for epsilon; 3: update epsilon according to likelihood (for 2 and 3, pre-specified epsilon was used as initial value) 
	int M[(n+1)], Mcan[(n+1)], Maccept[(n+1)];
	double P[(k+1)], Pcan[(k+1)]; 
	int Paccept[(k+1)];
	double ll[(n+1)][(k+1)];
	double llcan[(n+1)][(k+1)];
	double gridA[26][51];
	double gridB[26][51];
	double dataA1[(n+1)][(k+1)];
	double dataA2[(n+1)][(k+1)];
	double Strue[(n+1)][(k+1)];
	double Strue_can[(n+1)][(k+1)];
	int Strue_accept[(n+1)][(k+1)];
	int c_accept=0;
	double c_can;
	double q1=0.0, q2=0.0;
	for (i=1;i<=n;i++){
		M[i]= M0[i-1];
		Mcan[i]= M[i];
		Maccept[i]= 0;
		for (j=1;j<=k;j++){
			dataA1[i][j]= A1[(i-1)*k+j-1]; 
			dataA2[i][j]= A2[(i-1)*k+j-1]; 
			Strue[i][j]= dataA1[i][j]/(dataA1[i][j]+dataA2[i][j]);
			Strue_can[i][j]= Strue[i][j];
			Strue_accept[i][j]=0;
			if (i==1) {
				P[j]= P0[j-1];
				Pcan[j]= P[j];
				Paccept[j]= 0;
			}
		}
	}
	
	for (i=2;i<=25; i++){
		for (j=1; j<=50;j++){
			
			gridA[i][j]=A[(i-2)*50+j-1];
			gridB[i][j]=B[(i-2)*50+j-1];
			if (j==1) {
				gridA[i][0]= gridA[i][j];
				gridB[i][0]= gridB[i][j];
			}
		}
	}

	time_t t1, t2;
	t1 = time(NULL); // time 1

	////output////
	char var_file[1000];
	sprintf(var_file, "%s/%s", path[0], file_index[0]); 
	FILE *V0 = fopen(var_file, "w");
	
	////MCMC////
	
	//calculate likelihood of initial values
	for (i=1;i<=n;i++){
		for (j=1;j<=k;j++){
			ll[i][j]= logLike(M[i], P[j], dataA1[i][j], dataA2[i][j], Strue[i][j], gridA, gridB, c);
			llcan[i][j]= ll[i][j];
		}
	}
	
	for (i=1;i<=iter;i++){
		if (i%(iter/10)==0) Rprintf("Iter %d out of %d\n", i, iter);
		//update M
		for (j=1; j<=n; j++){
			prime = rbinom(1,0.5);
			if (prime==0) prime=-1;
			Mcan[j]= M[j] + prime;
			if ((Mcan[j]<=max_moi) && (Mcan[j]>0)){
				sumcan=0;
				sumori=0;
				for (x=1;x<=k;x++){
					llcan[j][x]=logLike(Mcan[j], P[x], dataA1[j][x], dataA2[j][x], Strue[j][x], gridA, gridB, c);
					sumcan+=llcan[j][x];
					sumori+=ll[j][x];
				}
				//accept
				if (log(runif(0.0,1.0)) < (sumcan-sumori)) {
					M[j]=Mcan[j];
					Maccept[j]++;
					for (x=1;x<=k;x++){
						ll[j][x]=llcan[j][x];
					}
				}
				//reject
				else {
					Mcan[j]=M[j];
					for (x=1;x<=k;x++){
						llcan[j][x]=ll[j][x];
					}
				}
			}
			else {
				Mcan[j]=M[j];
			}
		}
		//update PP
		for (j=1; j<=k; j++){
			Pcan[j]= rnorm(P[j],varP);
			if ((Pcan[j]<1) && (Pcan[j]>0)){
				sumcan=0;
				sumori=0;
				for (y=1;y<=n;y++){
					llcan[y][j]=logLike(M[y], Pcan[j], dataA1[y][j], dataA2[y][j], Strue[y][j], gridA, gridB, c);
					sumcan+=llcan[y][j];
					sumori+=ll[y][j];
				}
				//accept
				if (log(runif(0.0,1.0)) < (sumcan-sumori)) {
					P[j]=Pcan[j];
					Paccept[j]++;
					for (y=1;y<=n;y++){
						ll[y][j]=llcan[y][j];
					}
				}
				//reject
				else {
					Pcan[j]=P[j];
					for (y=1;y<=n;y++){
						llcan[y][j]=ll[y][j];
					}
				}
			}
			else {
				Pcan[j]=P[j];
			}
		}

		if (c>0){
			//update St
			for (j=1;j<=n;j++){
				for (x=1;x<=k;x++){
					Strue_can[j][x]= rnorm(Strue[j][x],varS);
					q1= dnorm(Strue[j][x],Strue_can[j][x],varS,0);
					q2= dnorm(Strue_can[j][x],Strue[j][x],varS,0);
					if (Strue_can[j][x]<0) {
						q2= pnorm(0, Strue[j][x], varS,1,0);
						Strue_can[j][x]=0;
						q1= dnorm(Strue[j][x],Strue_can[j][x],varS,0);
					}
					else if (Strue_can[j][x]>1) {
						q2= 1- pnorm(1, Strue[j][x], varS,1,0);
						Strue_can[j][x]=1;
						q1= dnorm(Strue[j][x],Strue_can[j][x],varS,0);			
					}
					if (Strue[j][x]==0){
						q1= pnorm(0, Strue_can[j][x], varS,1,0);
					}
					else if (Strue[j][x]==1){
						q1= 1- pnorm(1, Strue_can[j][x], varS,1,0);
					}
					
					llcan[j][x]=logLike(M[j], P[x], dataA1[j][x], dataA2[j][x], Strue_can[j][x], gridA, gridB, c);
					//accept
					if ((llcan[j][x]-ll[j][x]+log(q1)-log(q2))>log(runif(0.0,1.0))){
						Strue[j][x]= Strue_can[j][x];
						ll[j][x]= llcan[j][x];
						Strue_accept[j][x]++;
					}
					//reject
					else{
						Strue_can[j][x]= Strue[j][x];
						llcan[j][x]= ll[j][x];
					}
				} 
			}
		}
		
		//update c
		if (err_method==2){
			c = runif(0.0,upper_bound_c);
		}
		if (err_method==3){
			c_can= rnorm(c,varE);
			if (c_can>=0){
				//calculate likelihood
				sumcan=0;
				sumori=0;
				for (y=1;y<=n;y++){
					for (x=1;x<=k;x++){
						llcan[y][x]=logLike(M[y], P[x], dataA1[y][x], dataA2[y][x], Strue[y][x], gridA, gridB, c_can);
						sumcan+=llcan[y][x];
						sumori+=ll[y][x];
					}
				}
				//accept
				if (log(runif(0.0,1.0)) < (sumcan-sumori)) {
					c= c_can;
					c_accept++;
					for (y=1;y<=n;y++){
						for (x=1;x<=k;x++){
							ll[y][x]=llcan[y][x];
						}
					}
				}
				//reject
				else {
					c_can= c;
					for (y=1;y<=n;y++){
						for (x=1;x<=k;x++){
							llcan[y][x]=ll[y][x];
						}
					}
				}
			}
			else{
				//reject
				c_can= c;
			} //end update c
		}
		//print this iteration
		fprintf(V0,"%d", i);
		for (x=1;x<=n;x++) fprintf(V0,"\t%d",  M[x]);
		for (x=1;x<=k;x++) fprintf(V0,"\t%.6f", P[x]);
		if (err_method==3) fprintf(V0,"\t%.6f", c);
		fprintf(V0,"\n");		
	}
	
	fprintf(V0, "total_acceptance");
	for (x=1;x<=n;x++) fprintf(V0,"\t%d",  Maccept[x]);
	for (x=1;x<=k;x++) fprintf(V0,"\t%d", Paccept[x]);
	if (err_method==3) fprintf(V0,"\t%d", c_accept);
	fprintf(V0,"\n");
	
	t2 = time(NULL); // time 2
	Rprintf("Time = %.2f s\n", difftime(t2, t1));
	
	fclose(V0);
	PutRNGstate();
}

