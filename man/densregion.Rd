\name{densregion}
\alias{densregion}
\alias{densregion.default}
\title{ Density regions }
\description{
  A density region uses shading to represent the uncertainty surrounding a
  continuously-varying quantity, such as a survival curve or
  a forecast from a time series.  The darkness of the shading is
  proportional to the (posterior, predictive or fiducial) density.
  This function adds a density region to an existing plot. 
}
\usage{
densregion(x, \dots)
densregion.default(x, y, z, pointwise=FALSE, nlevels=100, 
                   colmax=par("fg"), scale=1, gamma=1, \dots)
}
\arguments{
  \item{x}{ Suppose the continuously-varying quantity varies over a
    space S.  \code{x} is a vector of the points in S at which the full
    posterior / predictive / fiducial distribution will be evaluated.  }
  \item{y}{ Vector of ordinates at which the density of the distribution
    will be evaluated for every point in \code{x}. 
  }
  \item{z}{ Matrix of densities on the grid defined by \code{x} and
    \code{y}, with rows corresponding to elements of \code{x} and
    columns corresponding to elements of \code{y}.  }
  \item{pointwise}{ If \code{TRUE} then the maximum density at each
    \code{x}
    is shaded with \code{colmax} (default black), and the shading
    intensity is proportional to the density within each \code{x}.
    
    If \code{FALSE} then the maximum density \emph{over all} \code{x} is shaded
    with \code{colmax}, and the shading is proportional to the density over
    all \code{x}. }
  \item{nlevels}{ Number of distinct shades to use to illustrate
    the varying densities.  The default of 100 should result in a plot with
    smoothly-varying shading.   }
  \item{colmax}{ Colour to shade the maximum density, either as a built-in R
    colour name (one of \code{\link{colors}()}) or an RGB hex
    value. Defaults to \code{par("fg")} which is normally     
    \code{"black"}, or \code{"#000000"}. }
  \item{scale}{ Proportion of \code{colmax} to shade the maximum density, for example \code{scale=0.5} with \code{colmax="black"} for a mid-grey colour. }
  \item{gamma}{Gamma correction to apply to the colour palette, see \code{\link{denstrip}}.}
  \item{\dots}{Further arguments passed to or from other methods.}
}
\details{
  The plot is shaded by interpolating the value of the density
  between grid points, using the algorithm described by Cleveland (1993)
  as implemented in the \code{\link{filled.contour}} function.

  With \pkg{lattice} graphics, similar plots can be implemented using
  the \code{\link{contourplot}} or \code{\link{levelplot}} functions.
}
\references{
  Jackson, C. H.  (2008) Displaying uncertainty with
  shading. \emph{The American Statistician}, 62(4):340-347.
  \url{http://www.mrc-bsu.cam.ac.uk/personal/chris/papers/denstrip.pdf}

  Cleveland, W. S. (1993) \emph{Visualizing Data}. Hobart Press, Summit,
  New Jersey. 
}
\author{ Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk> }
\seealso{ \code{\link{densregion.survfit}}, \code{\link{densregion.normal}}, \code{\link{denstrip}} }
\examples{
## Predictive uncertainty around a hypothetical regression line

x <- 1:10
nx <- length(x)
est <- seq(0, 1, length=nx)
lcl <- seq(-1, 0, length=nx)
ucl <- seq(1, 2, length=nx)
se <- (est - lcl)/qnorm(0.975)

y <- seq(-2, 2, length=100)
z <- matrix(nrow=nx, ncol=length(y))
for(i in 1:nx)
  z[i,] <- dnorm(y, est[i], se[i])
plot(x, type="n", ylim=c(-1.5, 2.5))
densregion(x, y, z, colmax="darkgreen")
lines(x, est)
lines(x, lcl, lty=2)
lines(x, ucl, lty=2)
\dontshow{
densregion(x, y, z, nlevels=20, colmax="darkgreen")
densregion(x, y, z, scale=0.5, colmax="darkmagenta")
}
## or automatically choose the y points to evaluate the density 

plot(x, type="n", ylim=c(-1.5, 2.5))
densregion.normal(x, est, se, ny=50, colmax="darkgreen")
lines(x, est)
lines(x, lcl, lty=2)
lines(x, ucl, lty=2)
}
\keyword{ aplot }
\keyword{ color }
