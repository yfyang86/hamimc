\name{MH.framework}
\alias{MH.framework}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
M-H framework 
}
\description{
M-H framework 
}
\usage{
MH.framework(dproposalfunction, rproposalfunction, dtarget, ITER = 10000, p = 50, x.initial = 1:p, buinin = 5000, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dproposalfunction}{
%%     ~~Describe \code{dproposalfunction} here~~
}
  \item{rproposalfunction}{
%%     ~~Describe \code{rproposalfunction} here~~
}
  \item{dtarget}{
%%     ~~Describe \code{dtarget} here~~
}
  \item{ITER}{
%%     ~~Describe \code{ITER} here~~
}
  \item{p}{
%%     ~~Describe \code{p} here~~
}
  \item{x.initial}{
%%     ~~Describe \code{x.initial} here~~
}
  \item{buinin}{
%%     ~~Describe \code{buinin} here~~
}
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
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
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (dproposalfunction, rproposalfunction, dtarget, ITER = 10000, 
    p = 50, x.initial = 1:p, buinin = 5000, ...) 
{
    X.flow = matrix(0, nrow = ITER, ncol = p)
    X.flow[1, ] = x.initial
    pb <- txtProgressBar(min = 0, max = ITER, style = 3)
    for (iter in 1:(ITER - 1)) {
        X.t = X.flow[iter, ]
        X.t1 <- rproposalfunction(X.t)
        if (log(runif(1)) < (dtarget(X.t1) + dproposalfunction(X.t, 
            X.t1) - dtarget(X.t) - dproposalfunction(X.t1, X.t))) {
            X.flow[iter + 1, ] = X.t1
        }
        else {
            X.flow[iter + 1, ] = X.t
        }
        setTxtProgressBar(pb, iter)
    }
    close(pb)
    return(X.flow)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~MHframe }

