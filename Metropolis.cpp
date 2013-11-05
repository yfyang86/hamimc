//
//  main.cpp
//  gibbs_sampler
//
//  Created by Yifan Yang on 11/4/13.
//  Copyright (c) 2013 Yifan Yang. All rights reserved.
//
#include <stddef.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <getopt.h>
#include <iostream>
#include <math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_matrix.h>

double y[]={28,-44,29,30,26,27,22,23,33,16,24,29,24,40,21,31,34,-2,25,19,24,28,37,32,20,25,25,36,36,21,28,26,32,28,26,30,36,29,30,22,36,27,26,28,29,23,31,32,24,27,27,27,32,25,28,27,26,24,32,29,28,33,39,25,16,23};
//#define dt(x,mu,sigma,scale,log)

double likelihood(double *param){
    double re=0.;
    for (int i=0;i<66;i++){
        re+= log(gsl_ran_tdist_pdf((y[i]-param[0])/sqrt(param[1]),5.))-log(param[1])/2.;
    }
    return (re);
}
/*
 muprior = dnorm(mu, mean=mu, sd=100,log = T)
 sigprior = dgamma(sigma2, shape = .1,rate=.1,log = T)
 */

double prior(double *param){
    double re1=0.,re2=0.;
    re1=gsl_ran_gaussian_pdf(param[0],100);
    re2=gsl_ran_gamma_pdf(param[1], .1, .1);
    return (log(re1*re2));
}

double posterior(double *param){
        return (likelihood(param) + prior(param));
}

inline double rnorm(gsl_rng *r,double mean,double sigma){
    return gsl_ran_gaussian(r, sigma)+mean;
}

void  proposalfunction(gsl_rng *r,double *param,double *re){
    re[0]=rnorm(r,param[0],1.);
    re[1]=rnorm(r,param[1],.7);//gsl_ran_chisq(r, 10.*sqrt(param[1]));
}

double runif(gsl_rng *r){
    return (gsl_rng_uniform(r));
}

void run_metropolis_MCMC(double *startvalue, int iterations,gsl_matrix * chain){
   // gsl_matrix * chain = gsl_matrix_alloc (iterations+1, 2);
    gsl_matrix_set(chain, 0, 0, startvalue[0]);
    gsl_matrix_set(chain, 0, 1, startvalue[1]);
    double param[2]={0.,0.};
    double proposaled[2]={0.,0.};
    double probab=0.;
    
    const gsl_rng_type * T;
    gsl_rng * r;
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);
    
    for (int i=0;i<iterations;i++){
        param[0]=gsl_matrix_get(chain, i, 0);
        param[1]=gsl_matrix_get(chain, i, 1);
        proposalfunction(r,param,proposaled);
        probab = exp(posterior(proposaled) - posterior(param));
        if (runif(r)<probab){
            gsl_matrix_set(chain, i+1, 0, proposaled[0]);
            gsl_matrix_set(chain, i+1, 1, proposaled[1]);
        }else{
            gsl_matrix_set(chain, i+1, 0, gsl_matrix_get(chain, i, 0));
            gsl_matrix_set(chain, i+1, 1, gsl_matrix_get(chain, i, 1));
        }
    }
    gsl_rng_free(r);
}


int my_gsl_matrix_fprintf(FILE *stream,gsl_matrix *m)
{
    size_t rows=m->size1;
    size_t cols=m->size2;
    size_t row,col,ml;
    int fill;
    char buf[100];
    gsl_vector *maxlen;
    
    maxlen=gsl_vector_alloc(cols);
    for (col=0;col<cols;++col) {
        ml=0;
        for (row=0;row<rows;++row) {
            sprintf(buf,"%g",gsl_matrix_get(m,row,col));
            if (strlen(buf)>ml)
                ml=strlen(buf);
        }
        gsl_vector_set(maxlen,col,ml);
    }
    
    for (row=0;row<rows;++row) {
        for (col=0;col<cols;++col) {
            sprintf(buf,"%g",gsl_matrix_get(m,row,col));
            fprintf(stream,"%s",buf);
            fill=gsl_vector_get(maxlen,col)+1-strlen(buf);
            while (--fill>=0)
                fprintf(stream," ");
        }
        fprintf(stream,"\n");
    }
    gsl_vector_free(maxlen);
    return 0;
}


int main(int argc, const char * argv[])
{
    unsetenv("VECLIB_MAXIMUM_THREADS");
    int itt=10000;
    if (argc==2){
        itt=atoi(argv[1]);
    }
    if (itt<1000){
        std::cerr<<"Min 1000, change to 10000"<<std::endl;
        itt=10000;
    }
    double startvalue[2]={25.,1.};
    FILE * f = fopen ("/Users/yifanyang/test.dat", "wb");
    //int dudu=10000;
    gsl_matrix * chain=gsl_matrix_alloc (itt+1,2);
    run_metropolis_MCMC(startvalue,itt,chain);
    my_gsl_matrix_fprintf(f, chain);
    fclose (f);
    gsl_matrix_free (chain);
    
    return 0;
}

