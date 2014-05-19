\name{dtrans}
\alias{dtrans}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Transverse comoving distance
%%  ~~function to do ... ~~
}
\description{Transverse comoving distance in Mpc between
    position (ra, dec) and reference point at position (ra0, dec0)
    with redshift z0.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
dtrans(ra, dec, z0, ra0, dec0, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{ra}{Right ascension in decimal degrees
%%     ~~Describe \code{ra} here~~
}
  \item{dec}{Declination in decimal degrees
%%     ~~Describe \code{dec} here~~
}
  \item{z0}{reference redshift
%%     ~~Describe \code{z0} here~~
}
  \item{ra0}{reference right ascension
%%     ~~Describe \code{ra0} here~~
}
  \item{dec0}{reference declination
%%     ~~Describe \code{dec0} here~~
}
  \item{\dots}{cosmological parameters passed to \code{\link{dcos}}
%%     ~~Describe \code{\dots} here~~
}
}
\details{
    Either the positions (ra, dec) or (z0, ra0, dec0) can be vectors.
    Cosmological parameters must be named if they are used.
    See \code{\link{dcos}}.
}
\value{
    The comoving transverse separation in Mpc.
}
\author{
  M.L. Peck
}


\examples{
    ## A compact group of galaxies from SDSS
    
    ra <- c(53.055043, 53.053015, 53.055655, 53.053242)
    dec <- c(0.26295811, 0.26203222, 0.26371558, 0.26455499)
    z <- c(0.0860, 0.08751, 0.08686, 0.08636)
    dtrans(ra, dec, z[1], ra[1], dec[1]) ##distance from object 1
    
}
\keyword{Basics}
\keyword{manip}
