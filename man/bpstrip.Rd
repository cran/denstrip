\name{bpstrip}
\alias{bpstrip}
\title{Box-percentile strips}
\description{
  Box-percentile strips give a compact illustration of a distribution. 
  The width of the strip is proportional to the probability of a more
  extreme point.  This function adds a box-percentile strip to an existing
  plot.
}
\usage{
bpstrip(x, prob, at, width, horiz=TRUE, scale=1, limits=c(-Inf, Inf), 
        col="gray", border=NULL, lwd=par("lwd"), lty=par("lty"), 
        ticks=NULL, tlen=1, twd=par("lwd"), tty=par("lty"))
}
\arguments{
  \item{x}{ Either the vector of points at which the probability is
    evaluated (if \code{prob} supplied), or a sample from the distribution
    (if \code{prob} not supplied). }
  \item{prob}{ Probability, or cumulative density, of the distribution
    at \code{x}.  If \code{prob} is not supplied, this is estimated from
    the sample \code{x} using \code{\link[stats]{ecdf}(x)}.}
  \item{at}{ Position of the centre of the strip on the y-axis (if
    \code{horiz=TRUE}) or the x-axis (if \code{horiz=FALSE}).  }
  \item{width}{ Thickness of the strip at its thickest point, which will
    be at the median.  Defaults to 1/20 of the axis range. }
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
}
\details{
  The box-percentile strip looks the same as the \emph{box-percentile plot}
  (Esty and Banfield, 2003) which is a generalisation of the boxplot for
  summarising data.  However, \code{\link{bpstrip}} is intended for illustrating
  distributions arising from parameter
  estimation or prediction.  Either the distribution is known
  analytically, or an arbitrarily large sample from the distribution is
  assumed to be available via a method such as MCMC or bootstrapping.
  
  The function \code{\link[Hmisc]{bpplot}} in the \pkg{Hmisc}
  package can be used to draw vertical box-percentile plots of observed
  data.  
}
\references{ Jackson, C. H.  (2008) Displaying uncertainty with
  shading. \emph{The American Statistician}, in press. 
  \url{http://www.mrc-bsu.cam.ac.uk/personal/chris/papers/denstrip.pdf}

  Esty, W. W. and Banfield, J. D. (2003) The box-percentile
  plot. \emph{Journal of Statistical Software} \bold{8}(17).  
}
\author{ Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk> }
\seealso{ \code{\link{vwstrip}}, \code{\link{cistrip}}, \code{\link{denstrip}} }
\examples{
x <- seq(-4, 4, length=1000)
prob <- pnorm(x)
plot(x, xlim=c(-5, 5), ylim=c(-5, 5), xlab="x", ylab="x", type="n")
bpstrip(x, prob, at=1, ticks=qnorm(c(0.25, 0.5, 0.75)))

## Terminate the strip at specific outer quantiles
bpstrip(x, prob, at=2, limits=qnorm(c(0.025, 0.975)))
bpstrip(x, prob, at=3, limits=qnorm(c(0.005, 0.995)))

## Compare with density strip
denstrip(x, dnorm(x), at=0)

## Estimate the density from a large sample 
x <- rnorm(10000)
bpstrip(x, at=4)
\dontshow{
x <- seq(-4, 4, length=1000)
bpstrip(x, prob, at=-1, col="purple", lwd=2, lty=2)
bpstrip(x, prob, at=-2, border=NA, ticks=qnorm(c(0.25,0.5,0.75)),
        twd=3,tty=2)
bpstrip(x, prob, at=-4, border=NA, horiz=FALSE)
bpstrip(x, prob, at=-3, border=NA, horiz=FALSE, width=0.2)
bpstrip(x, prob, at=2, border=NA, horiz=FALSE, scale=0.5)
}
}
\keyword{aplot}
