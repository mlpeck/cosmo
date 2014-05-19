## cosmological distances

dcos <- function(z, H0=70, Omega.m=0.27, Omega.l=1-Omega.m) {
  
  c <- 299792.458

  nz <- length(z)

  # Hubble distance (MPc)

  dH <- c/H0

  Omega.k <- 1-Omega.m-Omega.l

  E <- function(z) sqrt(Omega.m*(1+z)^3+Omega.k*(1+z)^2+Omega.l)

  fc <- function(z) 1/E(z)
  ft <- function(z) 1/((1+z)*E(z))

  # comoving distance

  dC <- numeric(nz)

  for (i in 1:nz) dC[i] <- dH*integrate(fc, lower=0, upper=z[i])$value

  # transverse comoving distance and 
  # comoving volume out to z

  if (Omega.k == 0) {
      dM <- dC 
      Vc <- (4*pi/3)*dM^3
  } else if (Omega.k > 0) {
      dM <- dH/sqrt(Omega.k)*sinh(sqrt(Omega.k)*dC/dH)
      Vc <- (2*pi*dH^3/Omega.k)*
         (dM/dH*sqrt(1+Omega.k*(dM/dH)^2) -
         1/sqrt(Omega.k)*asinh(sqrt(Omega.k)*dM/dH))
  } else {
    dM <- dH/sqrt(-Omega.k)*sin(sqrt(-Omega.k)*dC/dH)
    Vc <- (2*pi*dH^3/Omega.k)*
        (dM/dH*sqrt(1+Omega.k*(dM/dH)^2) -
        1/sqrt(-Omega.k)*asin(sqrt(-Omega.k)*dM/dH))
  }

  # angular diameter distance

  dA <- dM/(1+z)

  # luminosity distance

  dL <- (1+z)*dM

  # distance modulus

  dm <- 5*log10(dL)+25
  
  # comoving volume element per (steradian * delta z)
  
  dVc <- dH*(1+z)^2*dA^2*fc(z)

  # light travel distance (M light years)

  dT <- numeric(nz)

  for (i in 1:nz) dT[i] <- 3.261564*dH*integrate(ft, lower=0, upper=z[i])$value

  return(list(dC=dC, dM=dM, dA=dA, dL=dL, dm=dm, dT=dT, dVc=dVc, Vc=Vc))
}

## Flux to luminosity over solar bolometric luminosity
## note: solar bolometric luminosity is 3.839 e33 erg/sec.
## note: I've also seen 3.827 e 33. Which is it??
## 1 Mpc = 3.086 e24 cm
## sdss fluxes are in units of 1e-17 erg/cm^2/sec
## multiplying & dividing all the powers of 10 produces the net .01

lum.sol <- function(flux, z, ...) {
  dL <- dcos(z, ...)$dL
  lum <- 4*pi*0.01*(3.085678)^2/3.839*flux*dL^2
  lum[lum<0] <- NA
  return(lum)
}

loglum.ergs <- function(flux, z, ...) {
  dL <- dcos(z, ...)$dL
  lum <- 4*pi*(3.085678)^2*flux*dL^2
  lum[lum<0] <- NA
  return(31+log10(lum))
}


## Transverse comoving distance from (ra, dec) to (z0, ra0, dec0)

dtrans <- function(ra, dec, z0, ra0, dec0, ...) {
  dm <- dcos(z0, ...)$dM
  ra <- pi/180*ra
  dec <- pi/180*dec
  ra0 <- pi/180*ra0
  dec0 <- pi/180*dec0
  dx <- cos(ra)*cos(dec)-cos(ra0)*cos(dec0)
  dy <- sin(ra)*cos(dec)-sin(ra0)*cos(dec0)
  dz <- sin(dec)-sin(dec0)
  return(dm * sqrt(dx^2+dy^2+dz^2))
}

# velocity relative to a reference z = zbar

vrel <- function(z, zbar) 299792.458*(z-zbar)/(1+zbar)