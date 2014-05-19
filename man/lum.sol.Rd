\name{lum.sol}
\alias{lum.sol}
\title{Flux to luminosity in Sols
}
\description{Convert observed flux from an object at redshift z
    to luminosity over Sun's bolometric luminosity.
}
\usage{
lum.sol(flux, z, \dots)
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
    Luminosity in units of the solar bolometric luminosity. 
    Values < 0 return \code{NA}.
}
\author{
    M.L. Peck
}

\seealso{
    \code{\link{loglum.ergs}}
}
\keyword{Basics}
\keyword{manip}
