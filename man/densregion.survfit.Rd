\name{densregion.survfit}
\alias{densregion.survfit}
\title{Density regions for survival curves}
\description{
  Adds a density region to a survival plot.  The shading of the region
  has darkness proportional to the fiducial density of the point.  
  This distribution is assumed to be normal with standard deviation
  calculated using the lower confidence limit stored in the survival
  curve object.
}
\usage{
\method{densregion}{survfit}(x, ny=20, \dots)
}
\arguments{
  \item{x}{ Survival curve object, returned by
    \code{\link{survfit}}.  Confidence intervals must have been
    calculated, using \code{conf.type}.  }
  \item{ny}{ Minimum number of points to calculate the density at for
    each event time.  The default of 20 should be sufficient to obtain
    smooth-looking plots. }
  \item{\dots}{Further arguments passed to \code{\link{densregion.default}}.}
}
\details{
  The density is calculated at a grid of points, and interpolated using
  the method referred to in \code{\link{densregion}}.
}  
\note{
  This approach can only illustrate one survival curve per
  plot. Ideas would be welcome for how to define the shading for
  multiple curves, when the shaded areas representing the confidence
  regions overlap.  
}
\references{ Jackson, C. H.  (2008) Displaying uncertainty with
  shading. \emph{The American Statistician}, 62(4):340-347.
  \url{http://www.mrc-bsu.cam.ac.uk/personal/chris/papers/denstrip.pdf}
}
\author{ Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk> }
\seealso{ \code{\link{densregion}}, \code{\link{densregion.normal}}, \code{\link{denstrip}}}
\examples{
library(survival)
fit <- survfit(Surv(time, status), data=aml, conf.type="log-log")
plot(fit, col=0)
densregion(fit)
lines(fit, lwd=3, conf.int=FALSE, lty=1)
lines(fit, lwd=3, conf.int=TRUE, lty=2)
\dontshow{densregion(fit, ny=40)
densregion(fit, pointwise=TRUE)
densregion(fit, nlevels=10)
densregion(fit, colmax="purple")
densregion(fit, scale=0.6)}
## Wider CIs based on log survival
fit <- survfit(Surv(time, status), data=aml, conf.type="log")
plot(fit, col=0) 
densregion(fit) # Big variation in maximum density 
densregion(fit, pointwise=TRUE, colmax="maroon4")
par(new=TRUE)
plot(fit)

## Narrower CIs based on untransformed survival.
## Normal assumption probably unrealistic
fit <- survfit(Surv(time, status), data=aml, conf.type="plain")
plot(fit, col=0) 
densregion(fit, pointwise=TRUE, colmax="darkmagenta")
par(new=TRUE)
plot(fit)

}
\keyword{ aplot }
\keyword{ survival }
