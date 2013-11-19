
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]


double sum(NumericVector y){
double re=0.;
for (int i=0; i<y.size();i++) re+=y(i);
return re;
}

double sum(NumericVector x,NumericVector y){
double re=0.;
for (int i=0; i<y.size();i++) re+=y(i)*x(i);
return re;
}


List mcmctdist(SEXP inputy,SEXP iteration){
Environment stats("package:stats");
RNGScope scope;
double nu=5.;NumericVector yy(inputy);
NumericVector W(yy.size(),0.);
int N=yy.size();int i,j;int n=Rcpp::as<int>(iteration);NumericMatrix mat(n, 2);
double mu=20., sigma2=0.1;
for (i=0;i<n;i++){
  for (j=0;j<N;j++){
  	W(j)=R::rgamma((nu+1.)/2.,1./((yy(j)-mu)*(yy(j)-mu)/2+nu*sigma2/2.) );
		}
	  mu=R::rnorm(sum(W,yy)/(sum(W)+1.0/10000), 1./sqrt(sum(W)+1./10000));
    sigma2=R::rgamma(N*(nu+1.)/2.+0.1,1./(nu/2*sum(W)+0.1));
        mat(i,0) = mu;
        mat(i,1) = sigma2;
	}
List rest;rest["mat"]=mat;rest["mat2"]=mat;

    return rest;             // Return to R
}