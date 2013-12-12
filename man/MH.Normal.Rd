\name{MH.Normal}
\alias{MH.Normal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
H-H sampling for 50-variate Normal distribution.
}
\description{
H-H sampling for 50-variate Normal distribution.
}
\usage{
MH.Normal(ITER = 10000, p = 50, x.init = rep(0, 50), Sigma = matrix(rep(0.998, p^2), ncol = p) * outer(1:p, 1:p, "*") + 0.002 * diag((1:p)^2), chains = 4, buinin = 5000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ITER}{
Number of iterations
}
  \item{p}{
Number of dimensions
}
  \item{x.init}{
Initial value. Length=p.
}
  \item{Sigma}{
Var-Cov Matrix
}
  \item{chains}{
Number of chains
}
  \item{buinin}{
Buin-in steps. NO USE RIGHT now
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
	MH.Normal<-function(ITER=100);
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~mh }
