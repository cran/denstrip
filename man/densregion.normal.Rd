\name{densregion.normal}
\alias{densregion.normal}
\title{ Density regions based on normal distributions }
\description{
  Adds a density region to an existing plot of a normally-distributed quantity with
  continuously-varying mean and standard deviation, such as a
  time series forecast.  Automatically computes a reasonable set of 
  ordinates to evaluate the density at, which span the whole forecast
  space. 
}
\usage{
densregion.normal(x, mean, sd, ny=20, \dots)
}
\arguments{
  \item{x}{ Suppose the continuously-varying quantity varies over a
    space S. \code{x} is a vector of the points in S at which the
    posterior / predictive / fiducial distribution will be evaluated.  }
  \item{mean}{ Vector of normal means at each point in \code{x}. }
  \item{sd}{ Vector of standard deviations at each point in \code{x}. }
  \item{ny}{ Minimum number of points to calculate the density at for
    each \code{x}.  The density is calculated for at least \code{ny} equally
    spaced normal quantiles for each point.  The density is actually
    calculated at the union over \code{x} of all such points, for each \code{x}. 
  }
  \item{\dots}{Further arguments passed to \code{\link{densregion}}.}
}
\details{
  The plot is shaded by interpolating the value of the density
  between grid points, using the algorithm described by Cleveland (1993)
  as implemented in the \code{\link{filled.contour}} function.
}
\references{
  Jackson, C. H.  (2008) Displaying uncertainty with
  shading. \emph{The American Statistician}, 62(4):340-347.  \url{http://www.mrc-bsu.cam.ac.uk/personal/chris/papers/denstrip.pdf}

  Cleveland, W. S. (1993) \emph{Visualizing Data}. Hobart Press, Summit,
  New Jersey. 
}
\author{ Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk> }
\seealso{ \code{\link{densregion}}, \code{\link{densregion.survfit}}, \code{\link{denstrip}} }
\examples{
## Time series forecasting

(fit <- arima(USAccDeaths, order = c(0,1,1),
              seasonal = list(order=c(0,1,1))))
pred <- predict(fit, n.ahead = 36)
plot(USAccDeaths, xlim=c(1973, 1982), ylim=c(5000, 15000))

## Compute normal forecast densities automatically (slow)

densregion.normal(time(pred$pred), pred$pred, pred$se, 
                  pointwise=TRUE, colmax="darkgreen")
lines(pred$pred, lty=2)
lines(pred$pred + qnorm(0.975)*pred$se, lty=3)
lines(pred$pred - qnorm(0.975)*pred$se, lty=3)

## Compute forecast densities by hand (more efficient) 

nx <- length(pred$pred)
y <- seq(5000, 15000, by=100)
z <- matrix(nrow=nx, ncol=length(y))
for(i in 1:nx)
    z[i,] <- dnorm(y, pred$pred[i], pred$se[i])
plot(USAccDeaths, xlim=c(1973, 1982), ylim=c(5000, 15000))
densregion(time(pred$pred), y, z, colmax="darkgreen", pointwise=TRUE)
lines(pred$pred, lty=2)
lines(pred$pred + qnorm(0.975)*pred$se, lty=3)
lines(pred$pred - qnorm(0.975)*pred$se, lty=3)


densregion(time(pred$pred), y+2000, z, colmax="darkblue", pointwise=TRUE)
}
\keyword{ aplot }
\keyword{ color }
