\name{sectioned.density}
\alias{sectioned.density}
\alias{panel.sectioned.density}
\title{ Sectioned density plots }
\description{
  Sectioned density plots (Cohen and Cohen, 2006) use shading and
  occlusion to give a compact illustration of a distribution, such as
  the empirical distribution of data. 
}
\usage{
sectioned.density(x, dens, at, width, offset, ny,
                  method=c("kernel","frequency"), nx, horiz=TRUE,
                  up.left = TRUE, colmax, colmin="white", gamma=1,
                  lattice=FALSE, \dots)
panel.sectioned.density(...)
}
\arguments{
  \item{x}{ Either the vector of points at which the density is
    evaluated (if \code{dens} supplied), or a sample from the distribution
    (if \code{dens} not supplied).}
  \item{dens}{ Density at \code{points}. If \code{dens} is not supplied,
    the density of the distribution underlying \code{x} is estimated using the method
    specified in \code{method}. }
  \item{at}{ Position of the bottom of the plot on the y-axis (if
    \code{horiz=TRUE}) or position of the right of the plot on the
    x-axis (if \code{horiz=FALSE}) (required).  }
  \item{ny}{ Number of fixed-width intervals for categorising the density. }
  \item{width}{ Width of individual rectangles in the plot.  Defaults to
    the range of the axis divided by 20. }
  \item{offset}{ Offset for adjacent rectangles. Defaults to
    \code{width/3}. }
  \item{method}{Method of estimating the density of \code{x}, when
    \code{dens} is not supplied.
    
    If \code{"kernel"} (the default)
    then kernel density estimation is used, via \code{\link[stats]{density}(x,\dots)}. 

    If \code{"frequency"} then the density is estimated as the relative
    frequency in a series of bins, as in Cohen and Cohen (2006).  This
    method is controlled by the number of data bins \code{nx}.
  }
  \item{nx}{ Number of data bins for the \code{"frequency"} density
    estimation method. The default uses Sturges' formula (see \code{\link{nclass.Sturges}}, \code{\link{hist}}).}
  \item{horiz}{ If \code{horiz=TRUE}, then the plot is horizontal and
    points upwards.  If \code{horiz=FALSE} then the plot is vertical
    and points leftwards, as the illustrations in Cohen and Cohen
    (2006). }
  \item{up.left}{ If changed to \code{FALSE}, then horizontal plots
    point downwards and vertical plots point rightwards.}
  \item{colmax}{ Darkest colour, either as a built-in R colour name (one
    of \code{\link{colors}()}) or an RGB hex value.  Defaults to
    \code{\link{par}("fg")} or its \pkg{lattice} equivalent, which is normally     
    \code{"black"}, or \code{"#000000"}. }
  \item{colmin}{ Lightest colour, either as a built-in R colour name (one
    of \code{\link{colors}()}) or an RGB hex value.  Defaults to white. }
  \item{gamma}{Gamma correction to apply to the colour palette, see \code{\link{denstrip}}.}
  \item{lattice}{Set this to \code{TRUE} to make \code{\link{sectioned.density}}
    a \pkg{lattice} panel function instead of a base graphics function. \cr
    \code{panel.sectioned.density(x,...)} is equivalent to
    \code{sectioned.density(x, lattice=TRUE, ...)}.}
  \item{\dots}{Additional arguments supplied to \code{\link[stats]{density}(x,\dots)}, if
    \code{method="kernel"}. }
}
\references{ Cohen, D. J. and Cohen, J. The sectioned density
  plot. \emph{The American Statistician} (2006) \bold{60}(2):167--174 }
\author{ Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk> (R implementation)}
\examples{
## Fisher's iris data 
## Various settings to change the look of the plot
hist(iris$Sepal.Length, nclass=20, col="lightgray")
sectioned.density(iris$Sepal.Length, at=0.2)
sectioned.density(iris$Sepal.Length, at=5)
sectioned.density(iris$Sepal.Length, at=10, width=0.5)
hist(iris$Sepal.Length, nclass=20, col="lightgray")
sectioned.density(iris$Sepal.Length, at=7, width=0.5,
                  offset=0.1, colmax="darkmagenta")
sectioned.density(iris$Sepal.Length, at=9, width=0.5,
                  offset=0.1, ny=15, colmin="lemonchiffon")
## frequency method less smooth than kernel density 
sectioned.density(iris$Sepal.Length, at=12, width=0.5, offset=0.1,
                  method="frequency")
sectioned.density(iris$Sepal.Length, at=13.5, width=0.5, offset=0.1,
                  method="frequency", nx=20)

## Illustrate a known distribution
x <- seq(-4, 4, length=1000)
dens <- dnorm(x)
plot(x, xlim=c(-5, 5), ylim=c(-5, 5), xlab="x", ylab="x", type="n")
sectioned.density(x, dens, ny=8, at=0, width=0.3)
sectioned.density(x, dens, ny=16, at=2, width=0.1)
sectioned.density(x, dens, at=-3, horiz=FALSE)
sectioned.density(x, dens, at=4, width=0.3, horiz=FALSE)

}
\keyword{aplot}
