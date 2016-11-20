\name{ascale}
\alias{ascale}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Angular diameter scale in kpc/arcsec at redshift z
}
\description{
Calculates the angular diameter scale corresponding to the angular diameter distance
at redshift z.
}
\usage{
ascale(z, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{redshift
%%     ~~Describe \code{z} here~~
}
  \item{\dots}{cosmological parameters, if needed
%%     ~~Describe \code{\dots} here~~
}
}
\details{
The default cosmological parameters assumed are H0=70 km/sec/Mpc, Omega_m=0.27, 
Omega_lambda=1-Omega_m.
}
\value{
The angular scale in kpc/arcsec
}
\references{
    Hogg, D.W., 2000, "Distance measures in cosmology,"
    \href{http://arxiv.org/abs/astro-ph/9905116}{http://arxiv.org/abs/astro-ph/9905116}.
    
}
\author{
M.L. Peck
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{dcos}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{ascale(.05)
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (z, ...) 
{
    dA <- dcos(z, ...)$dA * 1000
    dA * tan(pi/(180 * 3600))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{math}% use one of  RShowDoc("KEYWORDS")
