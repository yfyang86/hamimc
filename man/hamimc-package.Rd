\name{hamimc-package}
\alias{hamimc-package}
\alias{hamimc}
\docType{package}
\title{
What the package does (short line)
}
\description{
More about what it does (maybe more than one line)
~~ A concise (1-5 lines) description of the package ~~
}
\details{
\tabular{ll}{
Package: \tab hamimc\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-11-18\cr
License: \tab What Licence is it under ?\cr
}
~~ An overview of how to use the package, including the most important ~~
~~ functions ~~
}
\author{
Who wrote it

Maintainer: Who to complain to <yourfault@somewhere.net>
}
\references{
~~ Literature or other references for background information ~~
}
~~ Optionally other standard keywords, one per line, from file KEYWORDS in ~~
~~ the R documentation directory ~~
\keyword{ package }
\seealso{
~~ Optional links to other man pages, e.g. ~~
~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
}
\examples{
## Not run
######################
## M-H tuned
######################
ITER=10000
x.initial=1:50
buinin=5000
dproposalfunction<-function(xt1,xt){# We could try other proposals
  xt<-as.vector(xt)
 Si=diag(rep(1,length(xt)))/20+xt%*%t(xt)/60
 mvtnorm::dmvnorm(xt1,mean=xt,sigma=Si,log=T)
}

rproposalfunction<-function(xt){
  xt<-as.vector(xt)
  Si=diag(rep(1,length(xt)))/20+xt%*%t(xt)/60
 mvtnorm::rmvnorm(n=1,mean=xt,sigma=Si)
}

#target function
A = matrix(rep(0.998,50^2),ncol=50)*outer(1:50,1:50,'*')
diag(A)=(1:50)^2

dtarget<-function(x,...){
 p=length(x)
 mvtnorm::dmvnorm(x,mean=rep(0,p),sigma=A,log=T)
}

MH.framework(dproposalfunction,rproposalfunction,dtarget,A=A) ->  X.flow

library(coda)

xyplot(mcmc(X.flow[,c(2,45)]))
par(mfrow=c(1,2))
densplot(mcmc(X.flow[-(1:buinin),c(2,45)]))
par(mfrow=c(1,1))

plot(X.flow[-(1:buinin),c(2,45)],type='l')
rejectionRate(mcmc(X.flow[-(1:buinin),c(2,45)]))
#     var1      var2 
# 0.9819964 0.9819964 

######################
##  Gibbs
######################

######################
##  HMC-plain
######################


######################
##  HMC
##  tune epsilon
######################

}
