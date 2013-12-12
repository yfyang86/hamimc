\name{hmcBiNorm.countour}
\alias{hmcBiNorm.countour}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cotour Plot
}
\description{
Contour Plot for Bi-Normal distribution
}
\usage{
hmcBiNorm.countour(X.flow, buinin = 500, bi = c(2, 45))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X.flow}{
Iteration X dimension matrix
}
  \item{buinin}{
Buinin steps
}
  \item{bi}{
Choose two dimensions to draw the contour plot. e.g. c(2,45)
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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Not run: 
X.flows <- HMC.cluster()
X.flow<-X.flows$valuechain[[1]];
hmcBiNorm.countour(X.flow)
## End(**Not run**)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~hmccontour }
