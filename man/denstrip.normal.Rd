\name{denstrip.normal}
\alias{denstrip.normal}
\title{Density strip for a normal or log-normal distribution}
\description{
  Draws a density strip for a normal or log-normal distribution
  with the given mean and standard deviation, based on computing the
  density at a large set of equally-spaced quantiles. 
}
\usage{
denstrip.normal(mean, sd, log=FALSE, nx=1000, \dots)
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
  \item{\dots}{Further arguments passed to \code{\link{denstrip}},
    for example, \code{at} to position the strip on the y-axis.}
}
\seealso{
  \code{\link{denstrip}}
}
\examples{
x <- seq(-4, 4, length=10000)
dens <- dnorm(x)
plot(x, xlim=c(-5, 5), ylim=c(-1, 2), xlab="x", ylab="", type="n", axes=FALSE)
axis(1)
denstrip(x, dens, at=0, width=0.3)
denstrip.normal(0, 1, at=1, width=0.3)

### log-normal distribution
sdlog <- 0.5
x <- rlnorm(10000, 0, sdlog)
plot(x, xlim=c(0, 5), ylim=c(-2, 4), xlab="x", , ylab="", type="n",
     axes=FALSE)
axis(1)
abline(v=1, lty=2, col="lightgray")
denstrip(x, at=0, ticks=exp(-sdlog^2), width=0.4) # tick at theoretical maximum density 
denstrip(x, at=1, bw=0.1, ticks=exp(-sdlog^2), width=0.4)
denstrip.normal(0, sdlog, log=TRUE, at=3, nx=1000,
                ticks=exp(-sdlog^2), width=0.4)
}
\author{ Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk> }
\keyword{aplot}
