\name{vwstrip}
\alias{vwstrip}
\title{Varying-width strips}
\description{
  Varying-width strips give a compact illustration of a distribution. 
   The width of the strip is proportional to the density.  This function 
   adds a varying-width strip to an exising plot.
  
}
\usage{
vwstrip(x, dens, at, width, horiz=TRUE, scale=1, limits=c(-Inf, Inf), 
        col="gray", border=NULL, lwd=par("lwd"), lty=par("lty"), 
        ticks=NULL, tlen=1, twd=par("lwd"), tty=par("lty"),\dots)
}
\arguments{
  \item{x}{ Either the vector of points at which the density is
    evaluated (if \code{dens} supplied), or a sample from the distribution
    (if \code{dens} not supplied). }
  \item{dens}{ Density at \code{x}. If \code{dens} is not supplied,
    the density of the sample \code{x} is estimated by kernel density
    estimation, using \code{\link[stats]{density}(x,\dots)}.}
  \item{at}{ Position of the centre of the strip on the y-axis (if
    \code{horiz=TRUE}) or the x-axis (if \code{horiz=FALSE}).  }
  \item{width}{ Thickness of the strip at the maximum density, that is, the length of its
    shorter dimension.  Defaults to 1/20 of the axis range. }
  \item{horiz}{ Draw the strip horizontally (\code{TRUE}) or vertically (\code{FALSE}).}
  \item{scale}{ Alternative way of specifying the thickness of the
    strip, as a proportion of \code{width}. }
  \item{limits}{Vector of minimum and maximum values, respectively, at which to terminate the strip.}
  \item{col}{Colour to shade the strip, either as a built-in R
    colour name (one of \code{\link{colors}()}) or an RGB hex
    value, e.g. black is \code{"#000000"}.}
  \item{border}{Colour of the border, see \code{\link{polygon}}. Use
    \code{border=NA} to show no border.  The default, 'NULL', means to
    use 'par("fg")'.}
  \item{lwd}{Line width of the border (defaults to \code{\link{par}("lwd")}).}
  \item{lty}{Line type of the border (defaults to \code{\link{par}("lty")}).}
  \item{ticks}{ Vector of \code{x}-positions on the strip to draw tick
    marks, or \code{NULL} for no ticks. }
  \item{tlen}{Length of the ticks, relative to the thickness of the strip.}
  \item{twd}{ Line width of these marks (defaults to \code{\link{par}("lwd")}). }
  \item{tty}{ Line type of these marks (defaults to \code{\link{par}("lty")}). }
  \item{\dots}{Additional arguments supplied to \code{\link[stats]{density}(x,\dots)}, if
    the density is being estimated. }
}
\details{
  Varying-width strips look like \emph{violin plots}.  The difference is that
  violin plots are intended to summarise data, while
  \code{\link{vwstrip}} is
  intended to illustrate a distribution arising from parameter
  estimation or prediction.  Either the distribution is known
  analytically, or an arbitrarily large sample from the distribution is
  assumed to be available via a method such as MCMC or bootstrapping.

  Illustrating outliers is important for summarising data, therefore
  violin plots terminate at the sample minimum and maximum and superimpose
  a box plot (which appears like the bridge of a violin, hence the name). 
  Varying-width strips, however, are used to illustrate known
  distributions which may have unbounded support.  Therefore it is
  important to think about where the strips should terminate (the
  \code{limits} argument).  For example, the end points may illustrate 
 a particular pair of extreme quantiles of the distribution.
  
  The function \code{\link[vioplot]{vioplot}} in the \pkg{vioplot}
  package and \code{\link[lattice]{panel.violin}} in the \pkg{lattice}
  package can be used to draw violin plots of observed data.
}
\references{ Jackson, C. H.  (2008) Displaying uncertainty with
  shading. \emph{The American Statistician}, in press. 
  \url{http://www.mrc-bsu.cam.ac.uk/personal/chris/papers/denstrip.pdf}

  Hintze, J.L. and Nelson, R.D. (1998) Violin plots: a box plot -
  density trace synergism. \emph{The American Statistician}
  \bold{52}(2),181--184. 
}
\author{ Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk> }
\seealso{ \code{\link{denstrip}}, \code{\link{bpstrip}}, \code{\link{cistrip}}.  }
\examples{
x <- seq(-4, 4, length=10000)
dens <- dnorm(x)
plot(x, xlim=c(-5, 5), ylim=c(-5, 5), xlab="x", ylab="x", type="n")
vwstrip(x, dens, at=1, ticks=qnorm(c(0.025, 0.25,0.5, 0.75, 0.975)))

## Terminate the strip at specific outer quantiles
vwstrip(x, dens, at=2, limits=qnorm(c(0.025, 0.975)))
vwstrip(x, dens, at=3, limits=qnorm(c(0.005, 0.995)))

## Compare with density strip
denstrip(x, dens, at=0)

## Estimate the density from a large sample 
x <- rnorm(10000)
vwstrip(x, at=4)
\dontshow{
x <- seq(-4, 4, length=10000)
vwstrip(x, dens, at=-1, col="purple", lwd=2, lty=2,
ticks=qnorm(c(0.25,0.5, 0.75)), twd=3, tty=2)
vwstrip(x, dens, at=-2, border=NA)
vwstrip(x, dens, at=-4, border=NA, horiz=FALSE)
vwstrip(x, dens, at=-3, border=NA, horiz=FALSE, width=1)
vwstrip(x, dens, at=2, border=NA, horiz=FALSE, scale=2)
}
}
\keyword{aplot}
