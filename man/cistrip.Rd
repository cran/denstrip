\name{cistrip}
\alias{cistrip}
\title{Line drawings of point and interval estimates}
\description{
  Adds one or more points and lines to a plot, representing point and
  interval estimates. 
}
\usage{
cistrip(x, at, d, horiz=TRUE, pch = 16, \dots)
}
\arguments{
  \item{x}{ Either a vector of three elements corresponding to point estimate,
    lower limit and upper limit of the interval estimate, respectively,
    or a matrix with three columns representing point estimates, lower
    and upper limits. 
  }
  \item{at}{ Position of the line on the y-axis (if
    \code{horiz=TRUE}) or the x-axis (if \code{horiz=FALSE}).  }
  \item{d}{Length of the serifs at each end of the line.  Defaults to
    1/60 of the axis range.}
  \item{horiz}{ Draw the line horizontally (\code{TRUE}) or vertically (\code{FALSE}).}  
  \item{pch}{Character to draw at the point estimate, see
    \code{\link{points}}.  By default this is a small solid circle, \code{pch=16}.}
  \item{\dots}{Further arguments passed to the \code{\link{points}} and
    \code{\link{segments}} functions.  For example \code{lty,lwd} to set the
    style and thickness of the line, or \code{cex} to expand the central
    point as in meta-analysis forest plots.}
}
\seealso{ \code{\link{denstrip}}, \code{\link{vwstrip}}, \code{\link{bpstrip}} }
\examples{
## One estimate 
x <- c(0.1, -2, 2)
plot(0, type="n", xlim=c(-5, 5), ylim=c(-5, 5), xlab="", ylab="")
abline(h=0, lty=2, col="lightgray")
abline(v=0, lty=2, col="lightgray")
cistrip(x, at=-0.1)
cistrip(x, at=0.2, lwd=3, d=0.1)
cistrip(x, at=-4, horiz=FALSE, lwd=3, d=0.2)

## Double / triple the area of the central point, as in forest plots
cistrip(x, at=2, d=0.2, pch=22, bg="black")
cistrip(x, at=2.5, d=0.2, pch=22, bg="black", cex=sqrt(2))
cistrip(x, at=3, d=0.2, pch=22, bg="black", cex=sqrt(3))

## Several estimates
x <- rbind(c(0.1, -2, 2), c(1, -1, 2.3),
           c(-0.2, -0.8, 0.4), c(-0.3, -1.2, 1.5))
plot(0, type="n", xlim=c(-5, 5), ylim=c(-5, 5), xlab="", ylab="")
cistrip(x, at=1:4)
abline(v=0, lty=2, col="lightgray")
cistrip(x, at=1:4, horiz=FALSE, lwd=3, d=0.2)
abline(h=0, lty=2, col="lightgray")
}
\author{ Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk> }
\keyword{aplot}
