\name{denstrip.legend}
\alias{denstrip.legend}
\alias{panel.denstrip.legend}
\title{ Add a legend to a density strip or shaded region }
\description{
  Add a legend to an existing plot with a density strip or shaded
  region, indicating the mapping of colours to densities.  
}
\usage{
denstrip.legend(x, y, width, len, colmax, gamma=1, horiz=FALSE, max, 
                nticks = 5, ticks, value.adj = 0, cex, main = "Density", lattice=FALSE)
panel.denstrip.legend(...)
}
\arguments{
  \item{x}{ Central x position of the legend. }
  \item{y}{ Central y position of the legend. }
  \item{width}{ Width of the legend strip, that is, the length of its
    shorter dimension. Defaults to 1/30 of the axis range.}
  \item{len}{ Length of the legend strip, that is, the length of its
    longer dimension. Defaults to 1/4 of the axis range.}
  \item{colmax}{ Colour at the maximum density, either as a built-in R
    colour name (one of \code{colors()}) or an RGB hex value.  Defaults
    to \code{par("fg")} or its \pkg{lattice} equivalent, which is
    \code{"black"} by default.}
  \item{gamma}{Gamma correction to apply to the colour palette, see \code{\link{denstrip}}.}
  \item{horiz}{ Legend strip drawn vertically (\code{FALSE}) or
    horizontally (\code{TRUE}).}
  \item{max}{ Maximum density on the legend, which is represented by \code{colmax}. }
  \item{nticks}{ Number of tick marks on the axis adjacent to the
    legend, if \code{ticks} not supplied. }
  \item{ticks}{Positions of numbered ticks on the axis adjacent to the
    legend. Defaults to \code{nticks} equally spaced ticks between 0 and
    the maximum density.}
  \item{value.adj}{ Extra adjustment for the axis labels to the right
    (if \code{horiz=FALSE}) or downwards (if \code{horiz=TRUE}). }
  \item{cex}{Text expansion.  Defaults to \code{\link{par}("cex") * 0.75} or
    \code{\link{trellis.par.get}("axis.text")$cex * 0.75}.}
  \item{main}{ Text to place above the legend. }
  \item{lattice}{Set this to \code{TRUE} to make \code{\link{denstrip.legend}}
    a \pkg{lattice} panel function instead of a base graphics function. \cr
    \code{panel.denstrip.legend(x,...)} is equivalent to
    \code{denstrip.legend(x, lattice=TRUE, ...)}.}
  \item{\dots}{Other arguments passed to \code{panel.denstrip.legend}.}
}
\author{ Christopher Jackson <chris.jackson@mrc-bsu.cam.ac.uk> }
\seealso{ \code{\link{denstrip}}, \code{\link{densregion}} }
\examples{
library(survival)
fit <- survfit(Surv(time, status) ~ 1, data=aml, conf.type="log-log")
plot(fit, col=0)
densregion(fit)
denstrip.legend(100, 0.8)
denstrip.legend(120, 0.8, width=3, len=0.4, value.adj=5)
denstrip.legend(40, 0.9, horiz=TRUE)
denstrip.legend(60, 0.7, horiz=TRUE, width=0.02, len=50, value.adj=0.04)
\dontshow{
denstrip.legend(40, 0.9, horiz=TRUE, colmax="darkmagenta")
denstrip.legend(60, 0.4, horiz=TRUE, colmax="darkmagenta", max=2)
}
}
\keyword{aplot}
\keyword{color}
