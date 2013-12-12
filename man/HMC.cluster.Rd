\name{HMC.cluster}
\alias{HMC.cluster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Run HMC multiple chains in parrellele.
}
\description{
Run HMC multiple chains in parrellele.
}
\usage{
HMC.cluster(nu.chains = 4, iterartion = 10000, startvalue = rep(0, 50))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nu.chains}{
Number of chains
}
  \item{iterartion}{
Number of iterations 
}
  \item{startvalue}{
Initial value.
}
}
\details{
Needs snow
}
\value{
mcmcchain is a mcmc.list in coda and valuechain store the value as a list. In each valuechain[[i]], there is a matrix store the multidimentional simulations.
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Yifan Yang
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#examples from HMC-clustered:
## Not run: 
system.time(X.flows <- HMC.cluster())
hmc.sum<-summary(X.flows$mcmcchain)
sum(hmc.sum$quantiles[,3]^2)
## End(**Not run**)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
