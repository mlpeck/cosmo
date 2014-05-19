\name{loglum.ergs}
\alias{loglum.ergs}
\title{Flux to log10 of luminosity in units ergs/sec
}
\description{Converts observed flux in units of 10^(-17) ergs/cm^2/sec
    to log10 of luminosity in ergs/sec
}
\usage{
loglum.ergs(flux, z, ...)
}
\arguments{
    \item{flux}{Observed flux in units of 10^(-17) erg/cm^2/sec}
    \item{z}{redshift}
    \item{\dots}{cosmological parameters passed to \code{\link{dcos}}}
}
\details{
    All values can be vectors.
    Cosmological parameters must be named if they are used.
    See \code{\link{dcos}}.
}
\value{
    log10 of the luminosity in ergs/sec. 
    Values < 0 return \code{NA}.
}
\author{
    M.L. Peck
}

\seealso{
    \code{\link{lum.sol}}
}
\keyword{Basics}
\keyword{manip}
