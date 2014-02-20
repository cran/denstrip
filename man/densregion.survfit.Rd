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
  In general, this approach can only illustrate one survival curve per
  plot.  Though if the graphics device supports transparency (e.g. PDF)
  multiple curves can be made to overlap smoothly - see the example
  below. 
}
\references{ Jackson, C. H.  (2008) Displaying uncertainty with
  shading. \emph{The American Statistician}, 62(4):340-347.
  \url{http://www.mrc-bsu.cam.ac.uk/personal/chris/papers/denstrip.pdf}
}
\author{ Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk> }
\seealso{ \code{\link{densregion}}, \code{\link{densregion.normal}}, \code{\link{denstrip}}}
\examples{
library(survival)
fit <- survfit(Surv(time, status) ~ 1, data=aml, conf.type="log-log")
plot(fit, col=0)
densregion(fit)
lines(fit, lwd=3, conf.int=FALSE, lty=1)
lines(fit, lwd=3, conf.int=TRUE, lty=2)
\dontshow{
plot(fit, col=0); densregion(fit, ny=40)
plot(fit, col=0); densregion(fit, pointwise=TRUE)
plot(fit, col=0); densregion(fit, nlevels=10)
plot(fit, col=0); densregion(fit, colmax="purple", new=FALSE)
plot(fit, col=0); densregion(fit, scale=0.6)
}
## Wider CIs based on log survival
fit <- survfit(Surv(time, status) ~ 1, data=aml, conf.type="log")
plot(fit, col=0) 
densregion(fit) # Big variation in maximum density 
plot(fit, col=0) 
densregion(fit, pointwise=TRUE, colmax="maroon4")
par(new=TRUE)
plot(fit)

## Narrower CIs based on untransformed survival.
## Normal assumption probably unrealistic
fit <- survfit(Surv(time, status) ~ 1, data=aml, conf.type="plain")
plot(fit, col=0) 
densregion(fit, pointwise=TRUE, colmax="darkmagenta")
par(new=TRUE)
plot(fit)

## Multiple survival curves on same axes
## Should overlap smoothly on devices that allow transparency
fit2 <- survfit(Surv(time, status) ~ x, data=aml, conf.type="log-log")
fit2x1 <- survfit(Surv(time, status) ~ 1, data=aml,
                  conf.type="log-log", subset=(x=="Maintained"))
fit2x0 <- survfit(Surv(time, status) ~ 1, data=aml,
                  conf.type="log-log", subset=(x=="Nonmaintained"))
plot(fit2, lwd=3, xlab="Weeks", ylab="Survival", xlim=c(0, 60),
     lty=1:2, col=c("red", "blue"), conf.int=TRUE, mark.time=TRUE)
densregion(fit2x1, colmax="red", gamma=2)
densregion(fit2x0, colmax="blue", gamma=2)

}
\keyword{ aplot }
\keyword{ survival }
