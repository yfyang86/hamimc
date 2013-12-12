\name{Gibbs.Norm50}
\alias{Gibbs.Norm50}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Gibbs sampling
}
\description{
Gibbs sampling for p dimension-Norm
}
\usage{
Gibbs.Norm50(ITER = 10000, p = 50, chains = 4)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ITER}{
Number of iterations
}
  \item{p}{
number of dimensions
}
  \item{chains}{
number of chains
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
\item{Gibbs.lst }{ returns a list contains chains of simulation result}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Zhiheng Xie, Yifan Yang
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
#examples from Gibbs:
## Not run: 
Gibbs.Norm50(1000)
## End(**Not run**)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
