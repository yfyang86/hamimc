\name{hmc.1}
\alias{hmc.1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
HMC main code
}
\description{
Hamiltonian Monte Carlo. This code is adapted from (Neal,2001). A simple tuning method is available.
}
\usage{
hmc.1(x, L = 1000, Tau = 0, epsilon = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
Initial value
}
  \item{L}{
Number of iterations
}
  \item{Tau}{
Steps of leapfrog. 
}
  \item{epsilon}{
Step-size of leapfrog.
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
Wang, Jian—Sheng, and R. H. Swendsen. "Nonuniversal critical dynamics in Monte Carlo simulations." Physical review letters (1987).
}
\author{

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
hmc.1(rep(0,50))
## End(**Not run**)	
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ HMC1 }

