\name{dcos}
\alias{dcos}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
    Cosmological distances and comoving volume
}
\description{
  Several measures of cosmological distances and comoving volume
  at redshift \code{z}.
}
\usage{
dcos(z, H0 = 70, Omega.m = 0.27, Omega.l = 1 - Omega.m)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{
    redshift
}
  \item{H0}{
    Hubble constant \eqn{H_0}{H0} in km/sec/Mpc.
}
  \item{Omega.m}{
    Dimensionless mass density parameter \eqn{\Omega_M}{Omega.m}
}
  \item{Omega.l}{
    Dimensionless cosmological constant parameter \eqn{\Omega_\Lambda}{Omega.l}
}
}
\details{\code{z} can be a vector.}
\value{
    A list with the following elements (units are Mpc unless otherwise noted):
    \item{dC}{Comoving distance}
    \item{dM}{Transverse comoving distance}
    \item{dA}{Angular diameter distance}
    \item{dL}{Luminosity distance}
    \item{dm}{Distance modulus (in magnitudes)}
    \item{dT}{Light travel time in Myrs.}
    \item{dVc}{Comoving volume element per (steradian * dz)}
    \item{Vc}{Comoving volume over the whole sky out to redshift z}
}
\references{
    Hogg, D.W., 2000, "Distance measures in cosmology,"
    \href{http://arxiv.org/abs/astro-ph/9905116}{http://arxiv.org/abs/astro-ph/9905116}.
}
\author{
    M.L. Peck
}

\seealso{
    \code{\link{dtrans}} for comoving transverse separation.
}
\examples{
    ## all distance measures at z=0.1
    
    dcos(0.1)
    
    ## reproduce figure 2 of Hogg
    
    z <- seq(0, 5, length=101)
    plot(z, dcos(z, Omega.m=1, Omega.l=0, H0=299792.458)$dA, type="l", ylim=c(0,0.5),
         xlab="redshift z", ylab="angular diameter distance dA/dH")
    points(z, dcos(z, Omega.m=0.05, Omega.l=0, H0=299792.458)$dA, type="l", lty=3)
    points(z, dcos(z, Omega.m=0.2, Omega.l=0.8, H0=299792.458)$dA, type="l", lty=2)
    
    ## reproduce figure 4 of Hogg
    
    z <- seq(0, 5, length=101)
    plot(z, dcos(z, Omega.m=1, Omega.l=0, H0=100)$dm, type="l", ylim=c(40,50),
         xlab="redshift z", ylab="distance modulus dm + 5 log h")
    points(z, dcos(z, Omega.m=0.05, Omega.l=0, H0=100)$dm, type="l", lty=3)
    points(z, dcos(z, Omega.m=0.2, Omega.l=0.8, H0=100)$dm, type="l", lty=2)
    
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{Basics}
\keyword{manip}
