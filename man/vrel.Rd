\name{vrel}
\alias{vrel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Velocity at redshift z relative to redshift zbar}
\description{
    The velocity at redshift z relative to a reference
    redshift zbar is simply \eqn{c (z-\bar{z})/(1+\bar{z})}{c * (z-zbar)/(1+zbar)}
}
\usage{
vrel(z, zbar)
}
\arguments{
  \item{z}{redshift}
  \item{zbar}{reference redshift}
}
\details{
    Either z or zbar (or both) can be vectors.
}
\value{Relative velocity in km/sec.}
\author{M.L. Peck}


\seealso{\code{\link{dtrans}} for transverse separation}
\examples{
    ## a compact group of galaxies from SDSS
    
    z <- c(0.0860, 0.08751, 0.08686, 0.08636)
    vrel(z, z[3])
}
\keyword{Basics}
\keyword{manip}
