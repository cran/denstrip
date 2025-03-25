\name{vwstrip.normal}
\alias{vwstrip.normal}
\alias{panel.vwstrip.normal}
\title{Varying width strip for a normal or log-normal distribution}
\description{
  Draws a varying width strip for a normal or log-normal distribution
  with the given mean and standard deviation, based on computing the
  density at a large set of equally-spaced quantiles. 
}
\usage{
vwstrip.normal(mean, sd, log=FALSE, nx=1000, \dots)
panel.vwstrip.normal(...)
}
\arguments{
  \item{mean}{Mean of the normal distribution.}
  \item{sd}{Standard deviation of the normal distribution.}
  \item{log}{If \code{TRUE} then the strip for a log-normal distribution,
    with mean and SD on the log scale \code{mean} and \code{sd},
    respectively, is plotted.  This may be useful for illustrating hazard ratios
    or odds ratios. 
  }
  \item{nx}{Number of points to evaluate the density at.}
  \item{\dots}{Further arguments passed to \code{\link{vwstrip}},
    for example, \code{at} to position the strip on the y-axis,
    or \code{lattice=TRUE} to use as a lattice panel function. \cr
  \code{panel.vwstrip.normal(x,\dots)} is equivalent to
  \code{vwstrip.normal(x, lattice=TRUE,\dots)}.}
}
\seealso{
  \code{\link{vwstrip}}
}
\examples{
x <- seq(-4, 4, length=10000)
dens <- dnorm(x)
plot(x, xlim=c(-5, 5), ylim=c(-1, 2), xlab="x", ylab="",
     type="n", axes=FALSE)
axis(1)
vwstrip(x, dens, at=0, width=0.4, limits=qnorm(c(0.005, 0.995)))
vwstrip.normal(0, 1, at=1, width=0.4, limits=qnorm(c(0.005, 0.995)))

### log-normal distribution
sdlog <- 0.5
x <- rlnorm(10000, 0, sdlog)
plot(x, xlim=c(0, 5), ylim=c(-1, 3), xlab="x", ylab="",
     type="n", axes=FALSE)
axis(1)
abline(v=1, lty=2, col="lightgray")
vwstrip(x, at=0, width=0.4, ticks=exp(-sdlog^2),
        limits=qlnorm(c(0.005,0.975),0,sdlog)) # tick at theoretical maximum density 
vwstrip(x, at=1, width=0.4, bw=0.1, ticks=exp(-sdlog^2),
        limits=qlnorm(c(0.005,0.975),0,sdlog))
vwstrip.normal(0, sdlog, log=TRUE, at=2.5, width=0.4, nx=1000,
               ticks=exp(-sdlog^2), limits=qlnorm(c(0.005,0.975),0,sdlog))
}
\author{ Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk> }
\keyword{aplot}
